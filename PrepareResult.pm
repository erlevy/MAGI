#!/usr/bin/perl
#04/09/2012
#Author: Petra Stepanowsky

package PrepareResult;

use warnings;
use strict;
use Statistics::R;
use File::Basename;
use File::Copy::Recursive qw(dircopy);
use File::Copy;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = qw(prepareresults);
@EXPORT_OK   = qw(prepareresults);
%EXPORT_TAGS = (DEFAULT  => [qw(&prepareresults)]);

##############################
# adds a directory path to files
# input: file names with no path, path of directory
# output: array of files names with full path
sub add_directory {
	my @files = @{$_[0]};
	my $dir = $_[1];
	for (my $i = 0; $i < scalar(@files); $i++) {
		$files[$i] = $dir.$files[$i];
	}
	return \@files;
}
##############################

##############################
# creates a list of reads mapped to given reference file
# input: mapped file, path of output directory, basename of sample (pattern), FASTA format of precursor file
# output: error message or 1
sub create_aln_file {
	#################
	# save input parameters
	my @files = @{$_[0]};
	my $output_dir = $_[1];
	my $precursors = $_[2];
	#################
	
	my %refs = ();
	
	if (-s $precursors) {
		# open reference file handler
		open (REF, "<", $precursors) or return "Cannot open file $precursors!\n";
		# save precursor ids and sequences
		while (<REF>) {
			my $line = $_;
			if ($line =~ /^>(.+)/) {
				my $ref = $1;
				my $seq = <REF>;
				$seq  =~ tr/[\r\n]//d; # remove \n and \r characters
				$refs{$ref} = $seq;
			}
		}
		# close reference file handler
		close REF;

		my %mapped = ();
		my %first_starts = ();
		my %last_ends = ();

		foreach my $file (@files) {
			# open input file handler
			open (IN, "<", $file) or return "Cannot open file $file!\n";
			my $basename = fileparse($file,  "_precursors.map");
			while(<IN>) {
				my $line = $_;
				if ($line =~ /^>(.+) count=(\d+) ref=(.+) strand=.+ start=(\d+) mismatches=\d+ cigar=(.*) align=(.+)/) {
					my $read = $1;
					my $count = $2;
					my $ref = $3;
					my $start = $4;
					my $cigar = $5;
					my $align = $6;
					my $seq = <IN>;
					$seq =~ tr/[\r\n]//d; # remove \n and \r characters
					
					# if alignment contains overhanging nts, subtract number of overhanging nts
					if ($cigar =~ /^(\d+)S/) {
						$start -= $1;
					}
					# only list reads without deletions or insertions
					if ($align !~ /\^/) {
						# calculate first possible start position of all reads per reference
						if (exists $first_starts{$ref}) {
							if ($first_starts{$ref} > $start) {
								$first_starts{$ref} = $start;
							}
						} else {
							$first_starts{$ref} = $start;
						}
						
						# calculate last possible end position of all reads per reference
						my $end = $start + length($seq);
					
						if (exists $last_ends{$ref}) {
							if ($last_ends{$ref} < $end) {
								$last_ends{$ref} = $end;
							}
						} else {
							$last_ends{$ref} = $end;
						}
						
						push (@{$mapped{$ref}{$basename}{$start}{$count}}, [$seq, $align]);
					} 
					# else skip reads with deletions
				}
			}
			# close input file handler
			close IN;
		}
		
		my $result = create_directory($output_dir."precursor_aln/"); return $result if ($result ne 1);
		
		foreach my $ref (sort keys %mapped) {
			open(OUT, ">", $output_dir."precursor_aln/".$ref.".aln") or return "Cannot open $ref.aln!\n";
			print OUT ">".$ref."\t".$refs{$ref}."\n";	
			foreach my $basename (sort keys %{$mapped{$ref}}) {				
				foreach my $start (sort {$a <=> $b} keys %{$mapped{$ref}{$basename}}) {
					foreach my $count (sort {$b <=> $a} keys %{$mapped{$ref}{$basename}{$start}}) {
						for (my $i = 0; $i < scalar(@{$mapped{$ref}{$basename}{$start}{$count}}); $i++) {
							print OUT $basename."\t".($start-1)."\t".$mapped{$ref}{$basename}{$start}{$count}[$i][0]."\t".$count."\n";
						}
					}
				}
			}
			close OUT;
		}
	}
	return 1;
}
##############################

##############################
# creates a correlation matrix of the microRNAs
# input: table formatted file, path of output directory
# output: 
sub create_correlation_result {
	#################
	# save input parameters
	my $norm_count_file = $_[0];
	my $output_dir = $_[1];
	#################
	
	if (-s $norm_count_file) {
		# create a communication bridge with R
		my $R = Statistics::R->new();

		# send variables to R
		$R->send(qq'file <- "$norm_count_file"');
		$R->send(qq'output_dir <- "$output_dir"');
		
		# prepare data
		$R->run(q'norm_count <- read.table(file, sep="\t", header=TRUE, row.names=1)');
		$R->run(q'da <- norm_count[,-ncol(norm_count)]');
		
		# calculate correlation matrix
		$R->run(q'cor.mat <- cor(da, method="spearman")');
		$R->run(q'cor.mat <- round(cor.mat, 4)');
		
		# write into a file
		$R->run(q'write.table(cor.mat, paste(output_dir, "norm_cor.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)');
		
		# stop the communication with R
		$R->stop();
	}
}
##############################

##############################
# creates a directory
# input: directory name/path
# output: error message or 1
sub create_directory {
	unless(-d $_[0]){
		mkdir $_[0] or return "Cannot create directory $_[0]!\n";
	}
	return 1;
}
##############################

##############################
# creates mapping report files per sample and read and for all samples
# input: read-count formatted files, mapped files with same basename as rc-files, path of output directory
# output: error message or 1
sub create_mapping_report {
	#################
	# save input parameters
	my @rc_files = @{$_[0]};
	my @mapped_files = @{$_[1]};
	my $output_dir = $_[2];
	#################
	
	my @report_files = ();
	
	my %sample_mapped = ();
	my %patterns = ();
	$patterns{"multiple"}++;
	
	my $multiple_total = 0;
	my $multiple_unique = 0;
	
	#################
	# for each read-count file save read information
	foreach my $file (@rc_files) {
		open (IN, "<", $file) or return "Cannot open file $file!\n";
		
		my $basename = fileparse($file, (".fa", ".rc", ".fasta"));
		my %reads_mapped = ();
		
		# save all reads in file in a hash
		while (<IN>) {
			my $line = $_;
			if ($line =~ /^>(.*) count=(\d+)/) {
				$reads_mapped{$1}{"count"} = $2;
				# total count
				$sample_mapped{$basename}{"reads"}[0] += $2;
				# unique count
				$sample_mapped{$basename}{"reads"}[1]++;
			}
		}
		close IN;
		
		# get all files with certain basename
		my @mapped = grep(/${basename}_/, @mapped_files);

		#################
		# for each basename corresponding mapped file save mapped read information
		foreach my $map (@mapped) {
			open (IN, "<", $map) or die "Cannot open file $map!\n";
			my $tmp = fileparse($map, ".map");
			my @splitted = split(/\_/, $tmp);
			my $pattern = $splitted[-1];
			$patterns{$pattern}++;
			# save reads mapped to certain database (pattern)
			while (<IN>) {
				my $line = $_;
				if ($line =~ /^>(.*) count=(\d+)/) {
					$reads_mapped{$1}{$pattern} = 1;
					# total count
					$sample_mapped{$basename}{$pattern}[0]+= $2;
					# unique count
					$sample_mapped{$basename}{$pattern}[1]++;
				}
			}
			close IN;
		}
		#################
	
		#################
		# write mapping report per sample if there are mapped reads
		if (scalar(keys %reads_mapped) > 0) {
			my $mappingreport_file = $output_dir.$basename."_mappingreport.txt"; 
			push (@report_files, $mappingreport_file);
			open (OUT, ">", $mappingreport_file) or return "Cannot create file $mappingreport_file!\n";
	
			# write header to file
			print OUT "reads\tcount";
			foreach my $pattern (sort keys %patterns) {
				print OUT "\t".$pattern;
			}
			print OUT "\n";
		
			# write read mapping information
			foreach my $key (sort {$a cmp $b} keys %reads_mapped) {
				print OUT $key."\t".$reads_mapped{$key}{"count"};
				
				my $flag_unique = 0;
				
				foreach my $pattern (sort keys %patterns) {
					if (exists $reads_mapped{$key}{$pattern}) {
						print OUT "\t".$reads_mapped{$key}{$pattern};
						if ($pattern ne "genome" && $pattern ne "mature") {
							$flag_unique += 1;
						}
					} else {
						print OUT "\t0";
					}
				}	
				print OUT "\n";
				
				if ($flag_unique > 1) {
					$sample_mapped{$basename}{"multiple"}[0] += $reads_mapped{$key}{"count"};
					$sample_mapped{$basename}{"multiple"}[1] += 1;
					foreach my $pattern (sort keys %patterns) {
						if ($pattern ne "genome" && exists $reads_mapped{$key}{$pattern}) {
							$sample_mapped{$basename}{$pattern}[0] -= $reads_mapped{$key}{"count"};
							$sample_mapped{$basename}{$pattern}[1] -= 1;
						}
					}
				}
			}
		
			close OUT;
		}
		#################
	}
	#################
	
	#################
	# write mapping report for all samples
	if (scalar(keys %sample_mapped) > 0) {
		open (TOTAL, ">", $output_dir."mapping_report_total.txt") or return "Cannot create file ${output_dir}mapping_report_total.txt!\n";
		open (UNIQUE, ">", $output_dir."mapping_report_unique.txt") or return "Cannot create file ${output_dir}mapping_report_unique.txt!\n";
		
		# write header to file
		print TOTAL "sample\treads";
		print UNIQUE "sample\treads";
		foreach my $pattern (sort keys %patterns) {
			print TOTAL "\t".$pattern;
			print UNIQUE "\t".$pattern;
		}
		print TOTAL "\n";	
		print UNIQUE "\n";			
		
		# write sample mapping information
		foreach my $key (sort {$a cmp $b} keys %sample_mapped) {	
			print TOTAL $key."\t".$sample_mapped{$key}{"reads"}[0];
			print UNIQUE $key."\t".$sample_mapped{$key}{"reads"}[1];
			
			foreach my $pattern (sort keys %patterns) {
				if (exists $sample_mapped{$key}{$pattern}) {
					if ($pattern eq "precursors" && exists $sample_mapped{$key}{"mature"}) {
						print TOTAL "\t".($sample_mapped{$key}{$pattern}[0] - $sample_mapped{$key}{"mature"}[0]);
						print UNIQUE "\t".($sample_mapped{$key}{$pattern}[1] - $sample_mapped{$key}{"mature"}[1]);
					} else {
						# total count
						print TOTAL "\t".$sample_mapped{$key}{$pattern}[0];
						# unique count
						print UNIQUE "\t".$sample_mapped{$key}{$pattern}[1];
					}
				} else {
					print TOTAL "\t0\t0";
					print UNIQUE "\t0\t0";
				}
			}	
			print TOTAL "\n";
			print UNIQUE "\n";
		}
		
		close TOTAL;
		close UNIQUE;
	}
	#################
	
	return 1;
}
##############################

##############################
# gets all files from a given directory with a certain pattern
# input: path of directory, pattern (regex)
# output: array of all found files
sub get_files_from_dir {
	#open directory handler
	opendir(DIR, $_[0]) or die "Cannot open directory $_[0]!\n";
	my @tmp = readdir(DIR);
	@tmp = sort(grep(/$_[1]/, @tmp));
	#close directory handler
	closedir(DIR);
	return add_directory(\@tmp, $_[0]);
}
##############################

##############################
# prepares result files and generates plots based on existing files
# input: read-count formatted files, path of output directory, FASTA format of precursor file, 
# output: error message or 1
sub prepareresults {
	#################
	# save input parameters
	# read count files
	my @rc_files = @{$_[0]};
	# an existing output directory
	my $dir = $_[1];
	# fasta file of precursor sequences
	my $file_precursors = $_[2];
	#################
	
	#################
	# specify file names
	my $dir_downsized_summary_total = $dir."downsized/downsize_summary_total.txt";
	my $dir_downsized_summary_unique = $dir."downsized/downsize_summary_unique.txt";
	#################
	
	#################
	# specify output directories
	my $dir_final = $dir."visualization/";
	my $dir_final_downsized = $dir_final."preprocessing/";
	my $dir_final_mapping = $dir_final."alignment/";
	my $dir_final_known = $dir_final."known/";
	my $dir_final_novel = $dir_final."novel/";
	my $dir_final_target = $dir_final."target/";
	#################
	
	#################
	# save basename of input files
	my @basenames = ();
	foreach my $file (@rc_files) {
		my $basename = fileparse($file, (".rc", ".fa", ".fq", ".fasta", ".fastq"));
		push (@basenames, $basename);
	}
	#################
	
	# create output directory
	my $res = create_directory($dir_final); return $res if ($res ne 1);
	
	#################
	# create downsized visualization files
	if (-d $dir."downsized/") {
		$res = create_directory($dir_final_downsized); return $res if ($res ne 1);
		
		if (-s $dir_downsized_summary_total && -s $dir_downsized_summary_unique) {
			$res = create_downsized_summary($dir_downsized_summary_total, $dir_final_downsized."downsize_summary_total.txt"); return $res if ($res ne 1);
			$res = create_downsized_summary($dir_downsized_summary_unique, $dir_final_downsized."downsize_summary_unique.txt"); return $res if ($res ne 1);
		}
	
		if (-d $dir."downsized/qual/") {
			dircopy($dir."downsized/qual/", $dir_final_downsized) or return "Cannot copy quality files!\n";
		}
		if (-d $dir."downsized/rc/") {
			dircopy($dir."downsized/rc/", $dir_final_downsized) or return "Cannot copy read-count files!\n";
		}
	}
	#################
	
	#################
	# collect mapped files (mapped to filter database)
	my @mapped_files = ();
	if (-d $dir."filtered/mapped") {
		push (@mapped_files, @{get_files_from_dir($dir."filtered/mapped/", '.map')});
	}
	#################
	
	#################
	if (-d $dir."known/") {
		my @precursors_mapped = ();
		if (-d $dir."known/mapped/") {
			push (@mapped_files, @{get_files_from_dir($dir."known/mapped/", '.map')});
			push (@precursors_mapped, @{get_files_from_dir($dir."known/mapped/", '_precursors.map')});
		}
		
		# create output directories
		$res = create_directory($dir_final_known); return $res if ($res ne 1);
		
	    # create alignment file for each precursor
		$res = create_aln_file(\@precursors_mapped, $dir_final_known, $file_precursors); return $res if ($res ne 1);
		
		################
		if (-d $dir."known/diff_expr/") {
			my $norm_count_file = $dir."known/diff_expr/norm_count.txt";
			
			$res = create_directory($dir_final_known."diff_expr/"); return $res if ($res ne 1);
			
			# create correlation plot of most abundant microRNAs
			create_correlation_result($norm_count_file, $dir_final_known."diff_expr/");
			
			################
			# copy DESeq results
			dircopy($dir."known/diff_expr/", $dir_final_known."diff_expr/") or return "Cannot copy differential expression files!\n";
			################
		}
		################
	}
	################
	
	################
	# create reports of mapping
	
	# create output directory
	$res = create_directory($dir_final_mapping); return $res if ($res ne 1);
	
	# create mapping report files
	$res = create_mapping_report(\@rc_files, \@mapped_files, $dir_final_mapping); return $res if ($res ne 1);
	################
	
	################
	if (-d $dir."combined/novel/") {
		$res = create_directory($dir_final_novel); return $res if ($res ne 1);
	
		my $novel_output_file = $dir_final_novel."precursors.str";
		create_novel_result($dir."combined/novel/precursors.str", $novel_output_file);
		my $novel_fa_file = $dir."combined/novel/novel.fa";
		if (-s $novel_fa_file) {
			copy($novel_fa_file, $dir_final_novel."novel.fa") or return "Cannot copy novel.fa!\n";
		}
		################
		if (-d $dir."novel/diff_expr/") {
			my $norm_count_file = $dir."novel/diff_expr/norm_count.txt";
			
			$res = create_directory($dir_final_novel."diff_expr/"); return $res if ($res ne 1);
			
			# create correlation plot of most abundant microRNAs
			create_correlation_result($norm_count_file, $dir_final_novel);
			
			################
			# copy DESeq results
			dircopy($dir."novel/diff_expr/", $dir_final_novel."diff_expr/") or return "Cannot copy differential expression files!\n";
			################
		}
		################
		
		################
		if (-d $dir."combined/novel/target/") {
			$res = create_directory($dir_final_target); return $res if ($res ne 1);
			
			my $target_input_file = $dir."combined/novel/target/novel_targets.txt";
			my $target_output_file = "novel_vs_genes.txt";
			$res = create_target_result($target_input_file, $target_output_file, $dir_final_target); return $res if ($res ne 1);
			copy($dir."combined/novel/target/novel_targets.txt", $dir_final_target."targets.txt");
		}
		################
	}	
	################
	
	return 1;
}
##############################

##############################
sub create_novel_result {
	my $file = $_[0];
	my $ofile = $_[1];
	
	open (IN, "<", $file) or return "Cannot open $file!\n";
	open (OUT, ">", $ofile) or return "Cannot create $ofile!\n";
	
	print OUT "id\tref\tstart\tstrand\tmature_start\tstar_start\tseq\tstruct\n";
	
	while (<IN>) {
		my $header = $_;
		$header =~ />(.*) ref=(.*) strand=(.*) mature=(.*) star=(.*)/;
		print OUT $1."\t".$2."\t".$3."\t".$4."\t".$5."\t";
		
		my $seq = <IN>;
		$seq =~ tr/[\r\n]//d; # remove \n and \r characters
		print OUT $seq."\t";
		
		my $struct = <IN>;
		$struct =~ tr/[\r\n]//d; # remove \n and \r characters
		print OUT $struct."\n";
	}
	
	close IN;
	close OUT;

	return 1;
}
##############################

##############################
sub create_downsized_summary {
	my $file = $_[0];
	my $ofile = $_[1];

	open (IN, "<", $file) or return "Cannot open $file!\n";
	open (OUT, ">", $ofile) or return "Cannot create $ofile!\n";
	
	my $first_line = 1;
	
	while (<IN>) {
		my $line = $_;
		my @split = split(/\t/, $line);
		if ($first_line) {
			print OUT $split[0]."\t".$split[1]."\t".$split[2]."\t".$split[4]."\t".$split[5]."\t".$split[6]."\t".$split[7]."\t".$split[8];
			$first_line = 0;
		} else {
			print OUT $split[0]."\t".$split[1]."\t".$split[2]."\t".$split[4]."\t".$split[5]."\t".$split[6]."\t".($split[7]-$split[5])."\t".$split[8];
		}
	}
	
	close IN;
	close OUT;
	
	return 1;
}
##############################

##############################
sub create_target_result {
	my $ifile = $_[0];
	my $ofile = $_[1];
	my $output_dir = $_[2];
	my %hash_miRNA = ();
	my $mirna_dir = $output_dir."miRNA/";
	my $gene_dir = $output_dir."gene/";
	
	if (-s $ifile) {
		open (IN, "<", $ifile) or return "Cannot open $ifile!\n";
		
		while (<IN>) {
			my $line = $_;
			
			if ($line =~ /\/\/hit_info\tquery_id=(.*)\treference_id=(.*)\tscore=(.*)\tenergy=(.*)\tquery_start=(\d+)\tquery_end=(\d+)\tref_start=(\d+)\tref_end=(\d+)\t.*\taln_mirna=(.*)\taln_map=(.*)\taln_utr=(.*)/) {
				my $query_id = $1;
				my $ref_id = $2;
				my $score = $3;
				my $energy = $4;
				my $q_start = $5;
				my $q_end = $6;
				my $r_start = $7;
				my $r_end = $8;
				my $aln_mirna = $9;
				my $aln_map = $10;
				my $aln_utr = $11;
				
				$ref_id =~ /(NM\_.+)\_(.+)/;
				my $transcript = $1;
				my $gene = $2;
				
				$hash_miRNA{$query_id}{$gene}{$transcript}{$score}{"energy"} = $energy;
				$hash_miRNA{$query_id}{$gene}{$transcript}{$score}{"q_start"} = $q_start;
				$hash_miRNA{$query_id}{$gene}{$transcript}{$score}{"q_end"} = $q_end;
				$hash_miRNA{$query_id}{$gene}{$transcript}{$score}{"r_start"} = $r_start;
				$hash_miRNA{$query_id}{$gene}{$transcript}{$score}{"r_end"} = $r_end;
				$hash_miRNA{$query_id}{$gene}{$transcript}{$score}{"aln_mirna"} = $aln_mirna;
				$hash_miRNA{$query_id}{$gene}{$transcript}{$score}{"aln_map"} = $aln_map;
				$hash_miRNA{$query_id}{$gene}{$transcript}{$score}{"aln_utr"} = $aln_utr;
			} 	
		}	
		close IN;
		
		my $res = create_directory($mirna_dir); return $res if ($res ne 1);
		$res = create_directory($gene_dir); return $res if ($res ne 1);
		
		open (SUM, ">", $output_dir.$ofile) or return "Cannot open $output_dir$ofile!\n";
		print SUM "miRNA\tgene\ttotal_score\n";
		foreach my $miRNA (keys %hash_miRNA) {
			open (MIRNA, ">", $mirna_dir.$miRNA.".txt") or return "Cannot create $mirna_dir$miRNA.txt!\n";
			foreach my $gene (keys %{$hash_miRNA{$miRNA}}) {
				open (GENE, ">>", $gene_dir.$gene.".txt") or return "Cannot create $gene_dir$gene.txt!\n";
				print SUM $miRNA."\t".$gene."\t";
				my $total_score = 0;
				foreach my $transcript (keys %{$hash_miRNA{$miRNA}{$gene}}) {
					foreach my $score (keys %{$hash_miRNA{$miRNA}{$gene}{$transcript}}) {
						print MIRNA $miRNA."\t".$gene."\t";
						print MIRNA $transcript."\t";
						print MIRNA $score."\t";
						print MIRNA $hash_miRNA{$miRNA}{$gene}{$transcript}{$score}{"energy"}."\t";
						print MIRNA $hash_miRNA{$miRNA}{$gene}{$transcript}{$score}{"q_start"}."\t";
						print MIRNA $hash_miRNA{$miRNA}{$gene}{$transcript}{$score}{"q_end"}."\t";
						print MIRNA $hash_miRNA{$miRNA}{$gene}{$transcript}{$score}{"r_start"}."\t";
						print MIRNA $hash_miRNA{$miRNA}{$gene}{$transcript}{$score}{"r_end"}."\t";
						print MIRNA $hash_miRNA{$miRNA}{$gene}{$transcript}{$score}{"aln_mirna"}."\t";
						print MIRNA $hash_miRNA{$miRNA}{$gene}{$transcript}{$score}{"aln_map"}."\t";
						print MIRNA $hash_miRNA{$miRNA}{$gene}{$transcript}{$score}{"aln_utr"}."\n";
						
						print GENE $gene."\t".$miRNA."\t";
						print GENE $transcript."\t";
						print GENE $score."\t";
						print GENE $hash_miRNA{$miRNA}{$gene}{$transcript}{$score}{"energy"}."\t";
						print GENE $hash_miRNA{$miRNA}{$gene}{$transcript}{$score}{"q_start"}."\t";
						print GENE $hash_miRNA{$miRNA}{$gene}{$transcript}{$score}{"q_end"}."\t";
						print GENE $hash_miRNA{$miRNA}{$gene}{$transcript}{$score}{"r_start"}."\t";
						print GENE $hash_miRNA{$miRNA}{$gene}{$transcript}{$score}{"r_end"}."\t";
						print GENE $hash_miRNA{$miRNA}{$gene}{$transcript}{$score}{"aln_mirna"}."\t";
						print GENE $hash_miRNA{$miRNA}{$gene}{$transcript}{$score}{"aln_map"}."\t";
						print GENE $hash_miRNA{$miRNA}{$gene}{$transcript}{$score}{"aln_utr"}."\n";
						
						$total_score += $score;
					}
				}
				print SUM $total_score."\n";
				close GENE;
			}
			close MIRNA;
		}
		
		close SUM;
	}
	
	return 1;
}
##############################

#loaded ok				
1;