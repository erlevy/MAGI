#!/usr/bin/perl
#03/29/2012
#Author: Petra Stepanowsky

package NovelPred;

use warnings;
use strict;

use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

use File::Basename;
use Clone qw(clone);
use Statistics::R;

use Aligner;

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(novelprediction);
%EXPORT_TAGS = (DEFAULT  => [qw(&novelprediction)]);

my $cnt = 1;
my %lTriplet = ("A" => {"010" => 0, "001" => 0, "011" => 0},
			   "C" => {"001" => 0, "110" => 0, "111" => 0},
			   "G" => {"100" => 0, "001" => 0, "110" => 0});
my %mTriplet = ("A" => {"000" => 0, "100" => 0, "001" => 0, "111" => 0},
			   "C" => {"000" => 0, "100" => 0, "110" => 0, "001" => 0},
			   "G" => {"000" => 0, "100" => 0, "010" => 0, "001" => 0, "110" => 0, "111" => 0},
			   "U" => {"000" => 0, "001" => 0});
my %rTriplet = ("A" => {"100" => 0, "010" => 0, "110" => 0},
			   "C" => {"100" => 0, "111" => 0},
			   "G" => {"001" => 0},
			   "U" => {"100" => 0, "001" => 0});

sub novelprediction {
	my $input_file = $_[0];
	my $output_dir = $_[1];
	my $genome_fasta_file = $_[2];
	my $genome_mismatches = $_[3];
	my $path_classifier = $_[4];
	my $type_classifier = $_[5];
	my $cutoff_pos = $_[6];
	my $read_count = $_[7];
	
	my $pp_fasta_file = $output_dir."potential_precursors.fa";
	my $pp_structure_file = $output_dir."potential_precursors.str";
	my $dataset_file = $output_dir."dataset_precursors.txt";
	my $probs_precursors_file = $output_dir."probs_precursors.txt";
	my $pp_bowtie2_index = $output_dir."potential_precursors";
	my $pos_bowtie2_index = $output_dir."positive_maturestar";
	my $input_file_basename = fileparse($input_file, ".map");
	
	my $res = create_directory($output_dir); return $res if ($res ne 1);
	$cnt = 1;
	$res = read_input_file($input_file); return $res if ($res =~ /Cannot open/);
	my %read_information = %{$res};
	$res = excise_potential_precursor($genome_fasta_file, \%read_information, $pp_fasta_file); return $res if ($res ne 1);
	$res = calculate_structure($pp_fasta_file, $pp_structure_file); return $res if ($res ne 1);
	
	#################
	# save secondary structure and mfe of potential precursors (pp) in hash
	my %pp_struct = ();
	my $id = "";
	open (STRUCT, "<", $pp_structure_file) or return "Cannot open $pp_structure_file!\n";
	while(<STRUCT>) {
		my $line = $_;
		$line =~ tr/\n\r//d;
		if ($line =~ /^>(.*) (ref=.*)/) {
			$id = $1; # identifier
			$pp_struct{$id}[0] = $2; # alignment information
		} elsif ($line =~ /^[ACGTUacgtu]+/) {
			$line =~ tr/Uu/Tt/; # sequence
			$pp_struct{$id}[1] = $line;
		} elsif ($line =~ /^([\.|\(|\)]+)\s\((.*)\)/) {
			$pp_struct{$id}[2] = $1; # structure
			$pp_struct{$id}[3] = $2; # minimum free energy
		}
	}	
	close STRUCT;
	#################
	
	#################
	# create bowtie2 index of potential precursors
	$res = create_bowtie2_index($pp_fasta_file, $pp_bowtie2_index); return $res if ($res ne 1);
	#################
	
	#################
	# align input file (reads) to potential precursors
	my @tmp_input = ($input_file);
	$res = Aligner::map(\@tmp_input, $pp_bowtie2_index, "pp", $output_dir, $genome_mismatches, "-a --norc"); return $res if ($res ne 1);
	#################
	
	my $mapped_file = $output_dir."mapped/".$input_file_basename."_pp.map";
	my %map_hash = ();
	
	open (MAP, "<", $mapped_file) or return "Cannot open $mapped_file!\n";
	while (<MAP>) {
		my $line = $_;
		$line =~ tr/\n\r//d;
		if ($line =~ /^>(.*) count=(\d+) ref=(.*) (strand=.*)/) {
			my $id = $1;
			my $count = $2;
			my $ref = $3;
			my $align = $4;
			my $seq = <MAP>;
			$seq =~ tr/\n\n//d;
			
			$map_hash{$ref}{$count}{$id} = [$align, $seq];
		}
	}
	close MAP;
	
    #open (FA, ">", $output_dir."precursors.fa") or return "Cannot create ${output_dir}precursors.fa!\n";
	
	my %pre_hash=();
	
	open (my $dataset, ">", $dataset_file) or return "Cannot create $dataset_file!\n";
	print_header($dataset);
	
	foreach my $ref_id (keys %map_hash) {
		my @read_counts = sort {$b <=> $a} keys %{$map_hash{$ref_id}};
		my $best_count = $read_counts[0];
		
		my @reads_highest_rc = keys %{$map_hash{$ref_id}{$best_count}};
		
		my $best_read_id = $reads_highest_rc[0];
		
		if (scalar(@reads_highest_rc) > 1) {
			my $best_align = $map_hash{$ref_id}{$best_count}{$best_read_id}[0];
			foreach my $read_id (@reads_highest_rc) {
				my $align = $map_hash{$ref_id}{$best_count}{$read_id}[0];
				
				if (is_better($align, $best_align)) {
					$best_read_id = $read_id;
					$best_align = $align;
				}
			}
		}

		$map_hash{$ref_id}{$best_count}{$best_read_id}[0] =~ /start=(\d+)/;
		my $start_mature = $1 - 1; # start of bowtie2 output is 1-based

		
		my $length_mature = length($map_hash{$ref_id}{$best_count}{$best_read_id}[1]);
		
		# pm -> potential mature
		my $struct_pm = substr($pp_struct{$ref_id}[2], $start_mature, $length_mature);
		
		if ($struct_pm !~ /([\(]+.*[\)]+)/ && $struct_pm !~ /([\)]+.*[\(]+)/) { # mature contains '(' and ')'
			my $struct_pp = $pp_struct{$ref_id}[2];
			my $seq_pp = $pp_struct{$ref_id}[1];
			my $found = 0;
			my $beg = -1;
			my $end = -1;
			my $length_star = 0;
			my @chars = split(//, $struct_pp);
			
			if ($struct_pm =~ /[\(]+/) { # mature is on 5' end
				my $pairs = ($struct_pm =~ tr/\(//);
				my $sum = $pairs;
				
				# star is on 3' end
				for (my $i = length($struct_pm) + $start_mature; $i < scalar(@chars) && !$found; ++$i) {
					if (substr($struct_pp, $i, 1) eq "(") {
						$sum += 1;
					} elsif (substr($struct_pp, $i, 1) eq ")") {
						$sum -= 1;
					}
					if ($sum - $pairs == 0) {
						$beg = $i + 1;
					}
					if ($sum == 0) {
						$end = $i;
						$found = 1;
					}
				}
				$length_star = $end-$beg+1;
				
			} elsif ($struct_pm =~ /[\)]+/) { # mature is on 3' end
				my $pairs = ($struct_pm =~ tr/\)//);
				my $sum = $pairs;
				
				# star is on 5' end
				for (my $i = $start_mature-1; $i >= 0 && !$found; --$i) {
					if (substr($struct_pp, $i, 1) eq ")") {
						$sum += 1;
					} elsif (substr($struct_pp, $i, 1) eq "(") {
						$sum -= 1;
					}
					if ($sum - $pairs == 0) {
						$end = $i;
					}
					if ($sum == 0) {
						$beg = $i;
						$found = 1;
					}
				}
				$length_star = $end-$beg;				
			}	

			# filter mature and star to get pre-miRNA
			my $mature_seq = substr($seq_pp, $start_mature, $length_mature);
			my $star_seq = substr($seq_pp, $beg, $length_star);
			
			if (length($mature_seq) <=32 && length($star_seq) <= 32 &&
			    abs(length($mature_seq) - length($star_seq)) <= 6) { # difference in length of mature and star only 6 nts
				my $start_pre = $start_mature;
				my $end_pre = $beg + $length_star - 1;
				if ($start_pre > $end_pre) {
					$start_pre = $beg;
					$end_pre = $start_mature + $length_mature - 1;
				}
				
				my $pre_seq = substr($seq_pp, $start_pre, $end_pre-$start_pre + 1);
				my $pre_struct = substr($struct_pp, $start_pre, $end_pre-$start_pre + 1);
				
				# single loop in pre-miRNA structure and length of pre-microRNA between 41 an 180 (according to real human pre-microRNAs
				if ($pre_struct !~ /\).*\(/ && (41 <= length($pre_seq) && length($pre_seq) <= 180)) { 
					$pre_struct =~ /(\().*\((\.*)\).*(\))\.*/;
					my @stemloop = ($-[2] + $start_pre, $+[2]-1 + $start_pre);
					
					if ($start_mature + $length_mature - 1 < $stemloop[0] || $start_mature > $stemloop[1]) { # mature not in loop part	
						$pre_hash{$ref_id}{"align"} = $pp_struct{$ref_id}[0];
						$pre_hash{$ref_id}{"best_read_id"} = $best_read_id;
						$pre_hash{$ref_id}{"mature_start"} = $start_mature;
						$pre_hash{$ref_id}{"mature_seq"} = $mature_seq;
						$pre_hash{$ref_id}{"star_start"} = $beg;
						$pre_hash{$ref_id}{"star_seq"} = $star_seq;
						$pre_hash{$ref_id}{"seq"} = $pre_seq;
						$pre_hash{$ref_id}{"struct"} = $pre_struct;
						
						calc_features($dataset, $pre_seq, $pre_struct, $pp_struct{$ref_id}[3], $ref_id);
					}
				}
			}
		} 
		#else {
			#print "mature doesn't have correct structure: $struct_pm\n";
		#}
	}
	#close FA;
	close $dataset;
	
	$res = run_classification($dataset_file, $path_classifier, $type_classifier, $probs_precursors_file); return $res if ($res ne 1);
	
	open (POS, "<", $probs_precursors_file) or return "Cannot open $probs_precursors_file!\n";
	open (MATURE, ">", $output_dir."novel_tmp.fa") or return "Cannot create ${output_dir}novel_tmp.fa!\n";
	
	my $cnt = 0;
	
	while(<POS>) {
		my $line = $_;
		if ($line =~ /(.*)\t(\d+\.\d+)/) {
			$pre_hash{$1}{"probability"} = $2;
			if ($2 > $cutoff_pos) {
				my $pId = $1;
				(my $id = $pId) =~ s/p//;
				
				print MATURE ">novel".$id."\n";
				print MATURE $pre_hash{$pId}{"mature_seq"}."\n";
				
				print MATURE ">novel".$id."*\n";
				print MATURE $pre_hash{$pId}{"star_seq"}."\n";
				
				$cnt +=1;
			}
		}
	}
	
	close MATURE;
	close POS;
	
	if ($cnt > 0) {
	
		#################
		# create bowtie2 index of novel mature and star sequences
		$res = create_bowtie2_index($output_dir."novel_tmp.fa", $pos_bowtie2_index); return $res if ($res ne 1);
		#################
		
		#################
		# align input file (reads) to positive classified precursors
		$res = Aligner::map(\@tmp_input, $pos_bowtie2_index, "posNovel", $output_dir, $genome_mismatches, "--norc"); return $res if ($res ne 1);
		#################
		
		my $novel_mapped_file = $output_dir."mapped/".$input_file_basename."_posNovel.map";
		
		open (IN, "<", $novel_mapped_file) or return "Cannot open $novel_mapped_file!\n";
		my %tmp = ();
		while (<IN>) {
			my $line = $_;
			$line =~ /count=(.*) ref=(.*) strand/;	
			$tmp{$2} += $1;
		}
		close IN;
		
		open (STR, ">", $output_dir."precursors.str") or return "Cannot create ${output_dir}precursors.str!\n";
		open (NOVEL, ">", $output_dir."novel.fa") or return "Cannot create ${output_dir}novel.fa!\n";
		
		foreach my $id (keys %tmp) {
			if ($tmp{$id} >= $read_count) {
			
				$id =~ s/novel//;
				$id =~ s/\*//;
				my $pId = "p".$id;
				
				my $align = $pre_hash{$pId}{"align"};
				$align =~ /.*(ref=.*) (strand=.*) (start=\d+).*/;
				print STR ">pre".$id." ".$1." ".$2." ".$3." mature=".$pre_hash{$pId}{"mature_start"}." star=".$pre_hash{$pId}{"star_start"}."\n";
				print STR $pre_hash{$pId}{"seq"}."\n";
				print STR $pre_hash{$pId}{"struct"}."\n";
				
				print NOVEL ">novel".$id."\n";
				print NOVEL $pre_hash{$pId}{"mature_seq"}."\n";
				
				print NOVEL ">novel".$id."*\n";
				print NOVEL $pre_hash{$pId}{"star_seq"}."\n";
			}
		}
		
		close STR;
		close NOVEL;
	}
	return 1;
}

##############################
#
# input:
# output: 
sub run_classification {
	my $dataset_file = $_[0];
	my $classifier_path = $_[1];
	my $classifier_type = $_[2];
	my $result_file = $_[3];

	# create a communication bridge with R
	my $R = Statistics::R->new();
	
	$R->send(qq'result_file <- "$result_file"');
	$R->send(qq'data_file <- "$dataset_file"');
	$R->send(qq'path_classifier <- "$classifier_path"');
	
	$R->run(q'data <- read.table(data_file, header=TRUE, sep="\t", row.names=1)');
	
	if ($classifier_type eq "LR") {
		$R->run(q'load(file=paste(path_classifier, "model_LR.RData", sep=""))'); # 'model'
		$R->run(q'if (nrow(data) > 0) {probs <- predict(model, newdata=data, type="response")} else { probs <- "" }');
	
	} elsif ($classifier_type eq "SVM") {
		$R->run(q'library(e1071)');
		$R->run(q'load(file=paste(path_classifier, "model_SVM.RData", sep=""))'); # 'model'
		$R->run(q'if (nrow(data) > 0) {probs <- attr(predict(model, newdata=data, probability=TRUE), "probabilities")[,1]} else { probs <- "" }');
	
	} elsif ($classifier_type eq "RF") {
		$R->run(q'library(randomForest)');
		$R->run(q'load(file=paste(path_classifier, "model_RF.RData", sep=""))'); # 'model'
		$R->run(q'if (nrow(data) > 0) {probs <- predict(model, newdata=data, type="prob")[,"pos"]} else { probs <- "" }');
		
	} else {
		return "Classifier method $classifier_type is not supported!\n";
	}
	
	$R->run(q'write.table(as.data.frame(probs), file=result_file, sep="\t", quote=FALSE, col.names=FALSE)');

	$R->stop();
	
	return 1;
}
##############################

##############################
#
# input:
# output: 
sub print_header {
	my $fh = $_[0];
	print $fh "id\tMFE_norm\tGC_ratio\tunpaired_paired_ratio\tnt_pairs";
		
	foreach my $nt (sort keys %lTriplet) {
		foreach my $t (sort {$a <=> $b} keys %{$lTriplet{$nt}}) {
			print $fh "\t".$nt.$t."_l";
		}
	}
	
	foreach my $nt (sort keys %mTriplet) {
		foreach my $t (sort {$a <=> $b} keys %{$mTriplet{$nt}}) {
			print $fh "\t".$nt.$t."_m";
		}
	}
	
	foreach my $nt (sort keys %rTriplet) {
		foreach my $t (sort {$a <=> $b} keys %{$rTriplet{$nt}}) {
			print $fh "\t".$nt.$t."_r";
		}
	}
}
##############################

##############################
#
# input: 
# output: 
sub calc_features {
	my $fh = $_[0];
	my $seq = $_[1];
	my $structure = $_[2];
	my $energy = $_[3];
	my $id = $_[4];
	
	$structure =~ /(\().*\((\.*)\).*(\))/;
	my $stembegin = $-[1];
	my $stemend = $-[3];
	my @stemloop = ($-[2], $+[2]-1);
	my $whole_stem_str = substr($structure, $stembegin, ($stemloop[0] - $stembegin)).substr($structure, $stemloop[1] + 1, ($stemend - $stemloop[1]));
		
	my %mTriplet = %{clone(\%mTriplet)};
	my %rTriplet = %{clone(\%rTriplet)};
	my %lTriplet = %{clone(\%lTriplet)};
	
	my $num_mtriplets = 0;
	my $num_ltriplets = 0;
	my $num_rtriplets = 0;

	for (my $i = 1; $i < length($seq) - 1; $i++) {
		my $structure_triplet = substr($structure, $i - 1, 3);
		$structure_triplet =~ s/\)/\(/g;
		$structure_triplet =~ s/\(/1/g;
		$structure_triplet =~ s/\./0/g;
		
		if (($i >= $stembegin) && ($i < $stemloop[0] || $i > $stemloop[1]) && ($i <= $stemend)) {
			#print substr($seq, $i, 1).$structure_triplet."\n";
			my $nt = substr($seq, $i, 1);
			$mTriplet{$nt}{$structure_triplet}++ if (exists $mTriplet{$nt}{$structure_triplet});
			$num_mtriplets++;
		}	

		if (($i > $stembegin) && ($i <= $stemloop[0] || $i > $stemloop[1] + 1) && ($i <= $stemend + 1)) {
			#print substr($seq, $i - 1, 1).$structure_triplet."\n";
			my $nt = substr($seq, $i - 1, 1);
			$lTriplet{$nt}{$structure_triplet}++ if (exists $lTriplet{$nt}{$structure_triplet});
			$num_ltriplets++;
		}	

		if (($i >= $stembegin - 1) && ($i < $stemloop[0] - 1 || $i >= $stemloop[1]) && ($i < $stemend)) {
			#print substr($seq, $i + 1, 1).$structure_triplet."\n";
			my $nt = substr($seq, $i + 1, 1);
			$rTriplet{$nt}{$structure_triplet}++ if (exists $rTriplet{$nt}{$structure_triplet});
			$num_rtriplets++;
		}					
	}			
		
	my $num_C = ($seq =~ tr/C//);
	my $num_G = ($seq =~ tr/G//);
	
	my $paired_nt = ($structure =~ tr/[\(\)]//);
	my $unpaired_nt = ($whole_stem_str =~ tr/\.//);
	my $nt_pairs = ($structure =~ tr/\(//);
	
	print $fh "\n".$id."\t";
	print $fh ($energy/length($seq))."\t";
	print $fh (($num_C + $num_G)/length($seq))."\t";
	print $fh ($unpaired_nt/$paired_nt)."\t";
	print $fh $nt_pairs;
	
	foreach my $nt (sort keys %lTriplet) {
		foreach my $triplet (sort {$a <=> $b} keys %{$lTriplet{$nt}}) {
			print $fh "\t".($lTriplet{$nt}{$triplet} / $num_ltriplets);
		}
	}
	foreach my $nt (sort keys %mTriplet) {
		foreach my $triplet (sort {$a <=> $b} keys %{$mTriplet{$nt}}) {
			print $fh "\t".($mTriplet{$nt}{$triplet} / $num_mtriplets);
		}
	}
	foreach my $nt (sort keys %rTriplet) {
		foreach my $triplet (sort {$a <=> $b} keys %{$rTriplet{$nt}}) {
			print $fh "\t".($rTriplet{$nt}{$triplet} / $num_rtriplets);
		}
	}
}
##############################

##############################
# 
# input: 
# output: 
sub is_better {
	my $a1 = $_[0];
	my $a2 = $_[1];

	$a1 =~ /mismatches=(\d+) cigar=(\d+)M/;
	my $a1_mis = $1;
	my $a1_matches = $2 - $a1_mis;
	
	$a2 =~ /mismatches=(\d+) cigar=(\d+)M/;
	my $a2_mis = $1;
	my $a2_matches = $2 - $a2_mis;
	
	if ($a1_matches > $a2_matches || # number of matches in a1 higher
	    ($a1_matches == $a2_matches && $a1_mis < $a2_mis)) { # number matches equal and number mismatches smaller in a1
		return 1;
	} 
	return 0;
}
##############################

##############################
# 
# input: 
# output: 
sub create_bowtie2_index {
	return (system("bowtie2-build $_[0] $_[1] 1>/dev/null 2>/dev/null") == 0);
}
##############################

##############################
# 
# input: 
# output: 
sub calculate_structure {
	return (system("RNAfold -noPS < $_[0] > $_[1]") == 0);
}
##############################

##############################
# 
# input: 
# output: 
sub excise_potential_precursor {
	my $genome_file = $_[0];
	my %read_info = %{$_[1]};
	my $o_file = $_[2];
	
	open (IN, "<", $genome_file) or return "Cannot open $genome_file!\n";
	open (my $fh, ">", $o_file) or return "Cannot open $o_file!\n";
	my $seq = "";
	my $id = "";
	my $description = "";
	while(<IN>) {
		my $line = $_;
		$line = clean_line($line);
		if ($line =~ /^>(\w+)(.*)/) {
			if ($id ne "" && $seq ne "") {
				my $res = excise_seq($id, $seq, \%read_info, $o_file, $fh); return $res if ($res ne 1);
			}
			$id = $1;
			$description = $2;	
			$seq = "";
		} else {
			$seq .= uc($line);
		}
	}
	my $res = excise_seq($id, $seq, \%read_info, $o_file, $fh); return $res if ($res ne 1);
	
	close IN;
	close $fh;
	
	return 1;
}
##############################

##############################
# 
# input: 
# output: 
sub excise_seq {
	my $id = $_[0];
	my $seq = $_[1];
	my %read_info = %{$_[2]};
	my $o_file = $_[3];
	my $fh = $_[4];
	
	foreach my $strand (sort keys %{$read_info{$id}}) {
		# most excised 3' position
		my $limit = 0;
		
		foreach my $start (sort {$a <=> $b} keys %{$read_info{$id}{$strand}}) {
			foreach my $end (sort {$a <=> $b} keys %{$read_info{$id}{$strand}{$start}}) {
				my $rc = $read_info{$id}{$strand}{$start}{$end}[0];
				my $rc_max = find_max_rc_downstream($start, $end, $read_info{$id}{$strand});
				
				if ($rc >= $rc_max && $start >= $limit) {
					my $start_precursor = $start - 70;
					my $end_precursor = $end + 20;
					
					my $pri_seq = substr($seq, $start_precursor, ($end_precursor - $start_precursor) + 1);
					if ($strand eq "-") {
						$pri_seq =~ tr/ACGTacgt/TGCAtgca/;
						$pri_seq = reverse($pri_seq);
					}
					
					print $fh ">p".$cnt." ref=".$id." strand=".$strand." start=".$start_precursor." end=".$end_precursor." read_id=".$read_info{$id}{$strand}{$start}{$end}[1]." read_count=".$rc." read_seq=".$read_info{$id}{$strand}{$start}{$end}[2]."\n";
					print $fh $pri_seq."\n";
					$cnt++;
					
					$start_precursor = $start - 20;
					$end_precursor = $end + 70;
					
					$pri_seq = substr($seq, $start_precursor, ($end_precursor - $start_precursor) + 1);
					if ($strand eq "-") {
						$pri_seq =~ tr/ACGTacgt/TGCAtgca/;
						$pri_seq = reverse($pri_seq);
					}
					
					print $fh ">p".$cnt." ref=".$id." strand=".$strand." start=".$start_precursor." end=".$end_precursor." read_id=".$read_info{$id}{$strand}{$start}{$end}[1]." read_count=".$rc." read_seq=".$read_info{$id}{$strand}{$start}{$end}[2]."\n";
					print $fh $pri_seq."\n";
					$cnt++;
					
					$limit = $end_precursor;
				}				
			}
		}
	}
	
	return 1;
}
##############################

##############################
# 
# input: 
# output: 
sub find_max_rc_downstream {
	my $start = $_[0];
	my $end = $_[1];
	my %hash_info = %{$_[2]};	
	my $max = 0;

	for (my $i = $start + 1; $i <= $end + 70; $i++) {
		if (exists $hash_info{$i}) {
			foreach my $j (sort {$a <=> $b} keys %{$hash_info{$i}}) {
				my $rc = $hash_info{$i}{$j}[0];
				if ($rc > $max) {
					$max = $rc;
				}
			}
		}
	}
	return $max;
}
##############################

##############################
# 
# input: 
# output: 
sub read_input_file {
	my $file = $_[0];
	my %hash_info = ();
	
	open (IN, "<", $file) or return "Cannot open $file!\n";
	
	my $read_id = "";
	my $expression = 0;
	my $ref = "";
	my $strand = "";
	my $start_pos = 0;
	
	while (<IN>) {
		my $line = $_;
		$line = clean_line($line);
		
		if ($line =~ /^>(.*)\scount=(\d+)\sref=(.*)\sstrand=(.?)\sstart=(\d+)/) {
			$read_id = $1;
			$expression = $2;
			$ref = $3;
			$strand = $4;
			$start_pos = $5;
		} else {
			my $seq = $line;
			my $end_pos = $start_pos + length($seq) - 1;
			$hash_info{$ref}{$strand}{$start_pos}{$end_pos} = [$expression, $read_id, $seq];
		}
	}
	
	close IN;
	
	return \%hash_info;
}
##############################
	
##############################
# Removes \n and \r from a string
# input: string
# output: string without \n and \r
sub clean_line {
	my $text = $_[0];
	$text =~ s/[\r\n]//g;
	return $text;
}
##############################

##############################
# Creates a directory
# input: directory name/path
# output: error message or 1
sub create_directory {
	unless(-d $_[0]){
		mkdir $_[0] or return "Cannot create directory: ($_[0])\n";
	}
	return 1;
}
##############################

#loaded ok				
1;