#!/usr/bin/perl
# v1.0
# 03/29/2012
# Author: Petra Stepanowsky

use strict;
use warnings;
use Getopt::Long;
use Config::Properties;
use File::Basename;
use File::Copy;
use Cwd 'abs_path', 'cwd';

use DownsizeFQ;
use Aligner;
use NovelPred;
use TargetPred;
use DiffExpr;
use PrepareResult;

my $usage = << "USAGE";
Description: 
Author: Petra Stepanowsky
Usage: perl .pl [options]
Options:
  --help,           --h              help
  --input,          --in     <file>  file in FASTQ format, more than one file separated by ','
  --dir,            --d      <dir>   directory containing FASTQ or read-count formatted files with extension '.fq','.fastq','.rc','.fasta' or '.fa'
  --outputdir,      --out    <dir>   output directory
  --config,         --c      <file>  Config properties file
Examples: perl .pl --input sample.fastq
          perl .pl --dir test/ 
USAGE

##### specify needed variables #####
my $help         =  0;

# path to bin and data directory
my $path_bin  = dirname(abs_path($0));
(my $path_data = $path_bin) =~ s/bin/data/;

# input variables
my $input_dir	 = "";
my $input 		 = "";
my @input_files  = ();

# path to config file 
my $config_file = $path_bin."/config.properties";
# current time and path of output directory
my $current_time = time();
my $output_dir   = cwd()."/pipeline_output".$current_time."/";

my %file_indexes = ();
my %file_aln = ();
#####

##### possible command options #####
GetOptions ('help|h'          => \$help,
			'dir|d=s'         => \$input_dir,
            'input|in=s'      => \$input,
			'outputdir|out=s' => \$output_dir,
			'config|c=s'      => \$config_file);
#####

##### print usage if the option 'help' is set or if no input directory or input file is specified #####
if ($help || ($input_dir eq "" && $input eq "")) {
	print $usage;
	exit;
}
#####

##### configuration variables #####
open (my $fh, "<", $config_file) or check_status("Cannot open file $config_file!\n");
my $config = Config::Properties->new();
$config->load($fh);
close $fh;
my $species = $config->getProperty("species");
# absolute path to specified species (data directory)
my $path_species = $path_data."/".$species."/";
#####

# read specified input files
get_input_files();

##############################
if (scalar(@input_files) > 0) {
	
	print "-----------------------------------\n";
	print "----- MiRAMAR analysis started -----\n";
	print "-----------------------------------\n";

	print "..... check database files: ";
	check_status(test_database_files());
	
	print "..... check input files:\n";
	my $file_format = test_input_files();
	check_status($file_format);
	
	print "..... create output directory: ";
	check_status(create_output_directory());
	
	print "..... write parameters: ";
	check_status(print_used_parameters());

##### downsize FASTQ files #####
	my @rc_files = ();
	if ($file_format eq "fastq") {
		print "..... downsize FASTQ files: ";
		
		my $res = DownsizeFQ::downsize(\@input_files, $output_dir, $config->getProperty("num_Ns_total"),
									   $config->getProperty("offset"), $config->getProperty("qual_type"),
									   $config->getProperty("qual_threshold"), $config->getProperty("adapter3"),
									   $config->getProperty("adapter_proportion"), $config->getProperty("num_mismatches"),
									   $config->getProperty("flag_count_mis_N"), $config->getProperty("min_length"),
									   $config->getProperty("max_length"), $config->getProperty("summary"),
									   $config->getProperty("qualfile"), $config->getProperty("trim"));
		
		return $res if ($res ne 1);
		
		@rc_files = @{get_files_from_dir($output_dir."downsized/rc/", '.rc')};
	} else {
		print "..... downsize FASTQ files: N/A\n";
		@rc_files = @input_files;
	}
#####
	
##### create combined folder and initialize file names for mapping and combined read file
	my $dir_combined = $output_dir."combined/";
	create_directory($dir_combined);
	
	my $read_mapping_file = $dir_combined."read_mapping.txt";
	my $combined_basename = "combined";
	my @combined_reads_file = ($dir_combined.$combined_basename.".rc");
#####	
	
##### combine all read-count files (if more than one file specified) to a unique set #####
	print "..... combine input files: ";
	if (scalar(@rc_files) == 1) {
		# copy input file to combined-directory (instead of combining)
		copy($rc_files[0], $combined_reads_file[0]) or check_status("Cannot copy $rc_files[0]!");
	} else {
		# combine reads and save mapping information
		combine_reads(\@rc_files, $read_mapping_file, $combined_reads_file[0]);
	}
	print "done\n";
#####
	
##### alignment to genome #####
	print "..... align to genome: ";
	my $res = Aligner::map(\@combined_reads_file, $file_indexes{"genome"}, "genome", $dir_combined, $config->getProperty("genome_mismatches"));
	check_status($res);
	# the genome mapped sequences are used for further alignments and process
	my @combined_genome_map_file = ($dir_combined."mapped/".$combined_basename."_genome.map");
	my $combined_genome_unmap_file = $dir_combined."unmapped/".$combined_basename."_genome.rc";
#####
	
##### filter reads #####
	print "..... filter reads: ";
	my @combined_filter_map_files = ();
	my @combined_filter_unmap_files = ();
	if ($config->getProperty("filter_reads") eq "yes") {
		# get all possible filter criteria
		opendir (DIR, $path_species."exclusioncriteria") or check_status("Cannot open species directory: $path_species\n");
		# get directory names of all filter criteria
		my @filter_dirs = readdir(DIR);
		# loop over each filter criteria
		print "\n";
		foreach my $dir_filter_name (@filter_dirs) {
			if ($dir_filter_name ne "." && $dir_filter_name ne "..") {
				print "     \t".$dir_filter_name.": ";
				# if there is a known index of the certain criteria available, align to filter sequences (information collected by testing databases)
				if (exists $file_indexes{$dir_filter_name}) {
					my $res = Aligner::map(\@combined_genome_map_file, $file_indexes{$dir_filter_name}, $dir_filter_name, $dir_combined, $config->getProperty($dir_filter_name."_mismatches"), "--norc");
					check_status($res);
					# save file name information to know filter mapping files
					push (@combined_filter_map_files, $dir_combined."mapped/".$combined_basename."_genome_".$dir_filter_name.".map");
					push (@combined_filter_unmap_files, $dir_combined."unmapped/".$combined_basename."_genome_".$dir_filter_name.".rc");
				} else {
					print "N/A\n";
				}
			}
		}
		closedir (DIR);
	} else {
		print "N/A\n";
	}
#####
	
##### detect known microRNAs #####	
	print "..... detect known microRNAs: ";
	my $combined_mature_map_file = "";
	my $combined_mature_unmap_file = "";
	my $combined_precursors_map_file = "";
	my $combined_precursors_unmap_file = "";
	my $combined_mature_survived_file = "";
	my @known_mapped_files = ();
	my @known_unmapped_files = ();
	if ($config->getProperty("detect_known") eq "yes") {
		# align to mature and mature star sequences
		print "\n     \talign to mature: ";
		my $res = Aligner::map(\@combined_genome_map_file, $file_indexes{"mature"}, "mature", $dir_combined, $config->getProperty("mature_mismatches"), "--local --norc");
		check_status($res);
		$combined_mature_map_file = $dir_combined."mapped/".$combined_basename."_genome_mature.map";
		push (@known_mapped_files, $combined_mature_map_file);
		$combined_mature_unmap_file = $dir_combined."unmapped/".$combined_basename."_genome_mature.rc";
		push (@known_unmapped_files, $combined_mature_unmap_file);
		
		# align to precursor sequences
		print "     \talign to precursors: ";
		$res = Aligner::map(\@combined_genome_map_file, $file_indexes{"precursors"}, "precursors", $dir_combined, $config->getProperty("precursors_mismatches"), "--norc");
		check_status($res);
		$combined_precursors_map_file = $dir_combined."mapped/".$combined_basename."_genome_precursors.map";
		push (@known_mapped_files, $combined_precursors_map_file);
		$combined_precursors_unmap_file = $dir_combined."unmapped/".$combined_basename."_genome_precursors.rc";
		push (@known_unmapped_files, $combined_precursors_unmap_file);
		
		# get only sequences that mapped to mature sequences and did not map to filter criteria
		print "     \tremove filtered reads: ";
		if ($config->getProperty("filter_reads") eq "yes") {
			my $res = get_survived_reads($dir_combined."survived/", $combined_mature_map_file, \@combined_filter_map_files, ".map");
			check_status($res);
			$combined_mature_survived_file = $dir_combined."survived/".$combined_basename."_genome_mature.map";
		} else {
			# if no filter criteria are specified, the survived reads are the mature mapped reads
			$combined_mature_survived_file = $combined_mature_map_file;
			print "N/A\n";
		}
	} else {
		print "N/A\n";
	}
#####
	
##### predict novel microRNAs #####
	print "..... predict novel microRNAs: \n";
	my $combined_genome_survived_file = "";
	if ($config->getProperty("predict_novel") eq "yes") {
		my @combined_exclude_files = (@combined_filter_map_files, $combined_mature_map_file, $combined_precursors_map_file);
		print "     \tremove reads mapped to other databases: ";
		my $res = get_survived_reads($dir_combined."survived/", $combined_genome_map_file[0], \@combined_exclude_files, ".map");
		check_status($res);
		print "     \tfind novel pre-miRNAs: ";
		$combined_genome_survived_file = $dir_combined."survived/".$combined_basename."_genome.map";
		
		$res = NovelPred::novelprediction($combined_genome_survived_file, $dir_combined."/novel/", $file_aln{"genome"}, $config->getProperty("genome_mismatches"), 
										  $path_bin."/", $config->getProperty("classifier"), $config->getProperty("cutoff_pos"), $config->getProperty("read_count"));
		check_status($res);
	} else {
		print "N/A\n";
	}
#####	

##### predict target sites for novel microRNAs #####
	my $novel_file = "";
	print "..... predict novel microRNA target sites: ";
	if ($config->getProperty("predict_novel") eq "yes" && $config->getProperty("predict_target") eq "yes") {
		create_directory($dir_combined."novel/target/");
		$novel_file = $dir_combined."novel/novel.fa";
		$res = TargetPred::targetprediction($novel_file, $file_aln{"mRNA"}, $config->getProperty("score"),
											$path_bin."/", $dir_combined."/novel/target/novel_targets.txt", 
											$config->getProperty("gpu"), $config->getProperty("rbs"), $config->getProperty("num_gpus"));
		check_status($res);
	} else {
		print "N/A\n";
	}
#####
		
##### split combined files #####
	print "..... split combined files: ";
	push (@combined_filter_map_files, $combined_genome_map_file[0]);
	push (@combined_filter_unmap_files, $combined_genome_unmap_file);
	if (scalar(@rc_files) == 1)	 {
		# copy combined files (mapped, unmapped and survived) and change basename to input file basename
		copy_files($rc_files[0], $combined_basename, \@combined_filter_map_files, \@combined_filter_unmap_files, \@known_mapped_files, 
					\@known_unmapped_files, $combined_mature_survived_file, $combined_genome_survived_file);
	} else {
		# split combined files to input files using the mapping file
		split_files($read_mapping_file, \@combined_filter_map_files, \@combined_filter_unmap_files, \@known_mapped_files, 
					\@known_unmapped_files, $combined_mature_survived_file, $combined_genome_survived_file);
	}
	print "done\n";
#####

	if ($config->getProperty("detect_known") eq "yes") {
		print "..... aggregate known microRNAs: ";
		check_status(aggregate_files(get_files_from_dir($output_dir."known/survived/", ".map"), $dir_combined."mapped/mature_list.txt", $output_dir."known/", "known"));
		my $aggregate_known_file = $output_dir."known/known_aggregate.bound";

##### diff. expr. for known microRNAs #####
		print "..... differential expression for known microRNAs:";
		if($config->getProperty("diff_expr_known") eq "yes") {
			my %group_names = ();
			
			foreach my $file (@rc_files) {
				my $basename = fileparse($file, (".rc", ".fasta", ".fa"));
				$basename =~ /(\w+?)\_/;
				$group_names{$1}++;
			}
			my @all_groups = sort keys %group_names;
			
			my @samples = ();
			my @groups = ();
			
			if (scalar(@all_groups) == 2) {
				foreach (my $i = 0; $i < scalar(@all_groups); $i++) {
					$samples[$i] = $group_names{$all_groups[$i]};
					$groups[$i] = $all_groups[$i];
				}
				print "\n     \t$groups[0] vs. $groups[1]: ";
				$res = 
				check_status(DiffExpr::diffexpression($aggregate_known_file, $output_dir."known/diff_expr/", \%group_names, \@groups, \@samples));
			
			} elsif (scalar(@all_groups) == 3) {
				# diff expression between group1 and group2
				$samples[0] = $group_names{$all_groups[0]};
				$groups[0] = $all_groups[0];
				$samples[1] = $group_names{$all_groups[1]};
				$groups[1] = $all_groups[1];
				print "\n     \t$groups[0] vs. $groups[1]: ";
				check_status(DiffExpr::diffexpression($aggregate_known_file, $output_dir."known/diff_expr/", \%group_names, \@groups, \@samples));
			
				# diff ecpression between group1 and group3
				$samples[1] = $group_names{$all_groups[2]};
				$groups[1] = $all_groups[2];
				print "     \t$groups[0] vs. $groups[1]: ";
				check_status(DiffExpr::diffexpression($aggregate_known_file, $output_dir."known/diff_expr/", \%group_names, \@groups, \@samples));
			
				# diff expression between group2 and group3
				$samples[0] = $group_names{$all_groups[1]};
				$groups[0] = $all_groups[1];
				print "     \t$groups[0] vs. $groups[1]: ";
				check_status(DiffExpr::diffexpression($aggregate_known_file, $output_dir."known/diff_expr/", \%group_names, \@groups, \@samples));
			} else {
				print "Differential expression only possible for 2 or 3 groups!\n";
			}
		} else {
			print "N/A\n";
		}
#####
	}

	if ($config->getProperty("predict_novel") eq "yes") {
		#print "..... aggregate novel microRNAs: ";
		#check_status(aggregate_files());
		
		#my $aggregate_novel_file = $output_dir."novel/novel_aggregate.bound";
		
##### diff. expr. for novel microRNAs #####
		print "..... differential expression for novel microRNAs: ";
		if($config->getProperty("diff_expr_novel") eq "yes") {
			#DiffExpr::diffexpression( , $output_dir."novel/diff_expr/", , "novel");
			print "done\n";
		} else {
			print "N/A\n";
		}
#####
	}

##### prepare and create results and plots #####
	print "..... create final report: ";
	$res = PrepareResult::prepareresults(\@rc_files, $output_dir, $file_aln{"precursors"});
	check_status($res);
#####

	print "------------------------------------\n";
	print "----- MiRAMAR analysis finished -----\n";
	print "------------------------------------\n";
} else {
	print "No input files found!\n";
}
##############################

##############################
# adds a directory path to files
# input: file names with no path, path of directory
# output: array of files names with full path			
sub add_directory {
	my @files = @{$_[0]};
	my $dir = $_[1];
	# for (my $i = 0; $i < scalar(@files); $i++) {
		# $files[$i] = $dir.$files[$i];
	# }
	foreach my $file (@files) {
		$file = $dir.$file;
	}
	return \@files;
}
##############################

##############################
# aggregates given files to one table formatted file with read count information of given reference ids
# input: array of map-formatted files, reference list, path of output directory, pattern
# output: error message or 1
sub aggregate_files {
	#################
	# save input parameters
	my @files = @{$_[0]};
	my $ref_file = $_[1];
	my $dir = $_[2];
	my $pattern = $_[3];
	#################
	
	my @samples = ();
	my %bound = ();
	
	#################
	# read and save reference information
	open(REF, "<", $ref_file) or return "Cannot open $ref_file!\n";
	my $flag = 1;
	while (<REF>) {
		my $line = $_;
		$line =~ tr/[\r\n]//d; # remove \n and \r characters
		if ($line =~ /^(.*)\tlength=\d+/) {
			my $id = $1;
			foreach my $file (@files) {
				my $basename = fileparse($file, ".map");
				$basename =~ /^([^_]+_[^_]+)/;
				$basename = $1;
				$bound{$id}{$basename} = 0;
				push (@samples, $basename) if ($flag);
			}
			$flag = 0;
		}
	}
	close REF;
	#################
	
	#################
	# add read count number per sample per reference id
	foreach my $file (@files) {
		my $basename = fileparse($file, ".map");
		$basename =~ /^([^_]+_[^_]+)/;
		$basename = $1;
	
		open(IN, "<", $file) or return "Cannot open $file!\n";
		while (<IN>) {
			my $line = $_;
			$line =~ tr/[\r\n]//d; # remove \n and \r characters
			if ($line =~ /^>.*\scount=(\d+)\sref=(.*)\sstrand=.+\sstart=\d+\smismatches=\d+\scigar=.*\salign=.*/) {
				if (exists $bound{$2}) {
					$bound{$2}{$basename} += $1;
				}
			}
		}
		close IN;
	}
	#################
	
	#################
	# print table formatted file
	open (OUT, ">", $dir.$pattern."_aggregate.bound") or return "Cannot create $dir${pattern}_aggregate.bound!\n";
	print OUT "miRNA";
	foreach my $sample (sort @samples) {
		print OUT "\t".$sample;
	}
	print OUT "\n";
	foreach my $miRNA (sort keys %bound) {
		print OUT $miRNA;
		foreach my $sample (sort keys %{$bound{$miRNA}}) {
			print OUT "\t".$bound{$miRNA}{$sample};
		}		
		print OUT "\n";
	}
	close OUT;
	#################
	
	return 1;
}
##############################

##############################
# checks the existence of a given database name
# input: path of directory, database name
# output: error message or 1
sub check_database {
	# get files from directory
	opendir (DIR, $_[1].$_[0]) or return "Cannot open directory $_[1]$_[0]!\n";
	my @dir_files = readdir(DIR);
	closedir(DIR);
	
	# select FASTA-formatted files (with extension .fa or .fasta)
	my @fasta_files = grep(/\.fasta$|\.fa$/, @dir_files);
	if (scalar(@fasta_files) != 1) {
		return "None or more than one possible $_[0] file found in $_[1]$_[0]!\n";
	} else {
		(my $basename = $fasta_files[0]) =~ s/\.fasta|\.fa//;
		$file_indexes{$_[0]} = $_[1].$_[0]."/".$basename;
		$file_aln{$_[0]} = $_[1].$_[0]."/".$fasta_files[0];
	}
	return 1;
}	
##############################	
	
##############################
# checks the bowtie2 indizes of given database pattern and directory
# input: path of directory, database name
# output: error message or 1
sub check_database_bowtie {
	# get files from directory
	opendir (DIR, $_[1].$_[0]) or return "Cannot open directory $_[1]$_[0]!\n";
	my @dir_files = readdir(DIR);
	closedir(DIR);
	
	# select FASTA-formatted files (with extension .fa or .fasta)
	my @fasta_files = grep(/\.fasta$|\.fa$/, @dir_files);
	# if it exists only one FASTA-formatted file, generate index files if they don't exist
	if (scalar(@fasta_files) == 1) {
		(my $basename = $fasta_files[0]) =~ s/\.fasta|\.fa//;
		$file_indexes{$_[0]} = $_[1].$_[0]."/".$basename;
		$file_aln{$_[0]} = $_[1].$_[0]."/".$fasta_files[0];
		my @bt_files = grep(/$basename.*\.bt2$/, @dir_files);
		unless (scalar(@bt_files) == 6) {
			#create bowtie2 index
			my $path = $_[1];
			my $dir = $_[0]."/";
			my $file = $fasta_files[0];
			my $res = system("bowtie2-build $path$dir$file $path$dir$basename 1>/dev/null 2>/dev/null");
			if ($res) {
				return "Cannot create bowtie2 index for $_[0] in $_[1]!\n";
			}
		}
		return 1;
	} else {
		return "None or more than one possible $_[0] file found in $_[1]$_[0]!\n";
	}
}
##############################

##############################
# checks the status of a given message
# input: result message of function call
# output: print "done" or call quit()
sub check_status {
	my $res = $_[0];
	if ($res eq 1) {
		print "done\n";
	} elsif ($res ne "fastq" && $res ne "rc" && $res !~ /Differential expression/) {
		print $res; 
		quit();
	}
}
##############################

##############################
# combines all reads from given files to a unique set of reads
# input: array of read-count formatted files, file name of mapping table, file name of combined read-count file
# output: NA
sub combine_reads {
	#################
	# save input parameters
	my @files = @{$_[0]};
	my $mapping_file = $_[1];
	my $combined_file = $_[2];
	#################
	
	my %hash = ();
	
	#################
	# foreach file saves information of header and sequence in hash
	foreach my $file (@files) {
		open (IN, "<", $file) or check_status("Cannot open file $file!\n");
		while (<IN>) {
			my $line = $_;
			$line =~ tr/[\r\n]//d; # remove \n and \r characters
			if ($line =~ /^>(.*)\scount=(\d+)/) {
				my $id = $1;
				my $rc = $2;
				my $seq = <IN>;
				$seq =~ tr/[\r\n]//d; # remove \n and \r characters
				$hash{$seq}{$file} = [$id, $rc];
			}
		}
		close IN;
	}
	#################

	#################
	# print header of mapping table
	open (OUT, ">", $mapping_file) or check_status("Cannot create file $mapping_file!\n");
	@files = sort (@files);
	foreach (my $i = 0; $i < scalar(@files); $i++) {
		my $basename = fileparse($files[$i], (".rc", ".fa", ".fasta"));
		print OUT "\t".$basename;
		print OUT "\t" if ($i + 1 < scalar(@files));
	}
	#################

	#################
	# print information to mapping table
	# print all reads to combined read-count file
	my $cnt = 1; # add number of processed reads as unique identifier
	open (OUTRC, ">", $combined_file) or check_status("Cannot create file $combined_file!\n");
	foreach my $seq (sort keys %hash) {
		my $sum_expr = 0;
		
		print OUT "\n".$seq;
		foreach my $file (sort @files) {
			if (exists $hash{$seq}{$file}) {
				print OUT "\t".$hash{$seq}{$file}[0]."\t".$hash{$seq}{$file}[1];
				$sum_expr += $hash{$seq}{$file}[1];
			} else { 
				print OUT "\t-\t0"; 
			}
		}
		
		print OUTRC ">read".$cnt." count=".$sum_expr."\n";
		print OUTRC $seq."\n";
		
		$cnt++;
	}
	close OUT;
	close OUTRC;
	#################
}
##############################

##############################
# copies the combined files to output directory and renames file name to sample name (if only one input file is specified)
# input: read-count file, basename of combined files, array of files mapped to filter databases, array of files unmapped to filter databases
# array of files mapped to known databases, array of files unmapped to known databases, survived known reads, survived filter reads
# output: NA
sub copy_files {
	#################
	# save input parameters
	my $rc_file = $_[0];
	my $combined_basename = $_[1];
	my @filter_files_map = @{$_[2]};
	my @filter_files_unmap = @{$_[3]};
	my @known_files_map = @{$_[4]};
	my @known_files_unmap = @{$_[5]};
	my $survived_known_file = $_[6];
	my $survived_filter_file = $_[7];
	#################
	
	my $dir_filtered = $output_dir."filtered/";
	my $dir_known = $output_dir."known/";
	
	#################
	#create directories
	if (scalar(@filter_files_map) > 0) {
		create_directory($dir_filtered);
		create_directory($dir_filtered."mapped/");
		create_directory($dir_filtered."unmapped/");
		create_directory($dir_filtered."survived/");
	}
	if (scalar(@known_files_map) > 0) {
		create_directory($dir_known);
		create_directory($dir_known."mapped/");
		create_directory($dir_known."unmapped/");
		create_directory($dir_known."survived/");
	}
	#################

	my $basename_rc_file = fileparse($rc_file, (".rc", ".fa", ".fasta"));
	#################
	# copy combined file to output directory and set name to sample name
	foreach my $filter_map_file (@filter_files_map) {
		(my $tmp_file_name = fileparse($filter_map_file)) =~ s/$combined_basename/$basename_rc_file/;
		copy($filter_map_file, $output_dir."filtered/mapped/".$tmp_file_name) or check_status("Cannot copy $filter_map_file to ${output_dir}filtered/mapped/$tmp_file_name!\n");
	}
	foreach my $filter_unmap_file (@filter_files_unmap) {
		(my $tmp_file_name = fileparse($filter_unmap_file)) =~ s/$combined_basename/$basename_rc_file/;
		copy($filter_unmap_file, $output_dir."filtered/unmapped/".$tmp_file_name) or check_status("Cannot copy $filter_unmap_file to ${output_dir}filtered/unmapped/$tmp_file_name!\n");
	}
	foreach my $known_map_file (@known_files_map) {
		(my $tmp_file_name = fileparse($known_map_file)) =~ s/$combined_basename/$basename_rc_file/;
		copy($known_map_file, $output_dir."known/mapped/".$tmp_file_name) or check_status("Cannot copy $known_map_file to ${output_dir}known/mapped/$tmp_file_name!\n");
	}
	foreach my $known_unmap_file (@known_files_unmap) {
		(my $tmp_file_name = fileparse($known_unmap_file)) =~ s/$combined_basename/$basename_rc_file/;
		copy($known_unmap_file, $output_dir."known/unmapped/".$tmp_file_name) or check_status("Cannot copy $known_unmap_file to ${output_dir}knoen/unmapped/$tmp_file_name!\n");
	}
	if ($survived_known_file ne "") {
		(my $tmp_file_name = fileparse($survived_known_file)) =~ s/$combined_basename/$basename_rc_file/;
		copy($survived_known_file, $output_dir."known/survived/".$tmp_file_name) or check_status("Cannot copy $survived_known_file to ${output_dir}known/survived/$tmp_file_name!\n");
	}
	if ($survived_filter_file ne "") {
		(my $tmp_file_name = fileparse($survived_filter_file)) =~ s/$combined_basename/$basename_rc_file/;
		copy($survived_filter_file, $output_dir."filtered/survived/".$tmp_file_name) or check_status("Cannot copy $survived_filter_file to ${output_dir}filtered/survived/$tmp_file_name!\n");
	}
	#################
}
##############################

##############################
# creates a directory
# input: directory name/path
# output: error message or 1		
sub create_directory {
	unless(-d $_[0]) {
		mkdir $_[0] or return "Cannot create directory $_[0]!\n"
	}
	return 1;
}	
##############################

##############################
# creates the output directory
# input: NA
# output: NA	
sub create_output_directory {
	return create_directory($output_dir);
}
##############################

##############################
# gets all files from a given directory with a certain pattern
# input: path of directory, pattern (regex)
# output: array of all found files
sub get_files_from_dir {
	my @tmp = ();
	if (-d $_[0]) {
		opendir(DIR, $_[0]) or check_status("Cannot open directory $_[0]!\n");
		@tmp = readdir(DIR);
		@tmp = sort(grep(/$_[1]/, @tmp));
		#close directory handler
		closedir(DIR);
		return add_directory(\@tmp, $_[0]);
	} 
	return \@tmp;
}
##############################

##############################
# reads the input files from a given directory (file extension .fastq, .fq or .fasta, .fa)
# or splits the given input string of input files
# input: NA
# output: NA
sub get_input_files {
	if ($input_dir ne "") {
		# open directory handler
		opendir(DIR, $input_dir) or check_status("Cannot open directory $input_dir!\n");
		my @dir_files = readdir(DIR);
		# close directory handler
		closedir(DIR);
		# try to get fastq files
		@input_files = sort(grep(/\.fastq$|\.fq$/, @dir_files));
		# if no fastq files found, try to get fasta files
		if (scalar (@input_files) <= 0) {
			@input_files = sort(grep(/\.fasta$|\.fa$|\.rc$/, @dir_files));
		}	
		@input_files = @{add_directory(\@input_files, $input_dir)};
	} elsif ($input ne "") {
		# split input file sequence
		@input_files = sort(split(/,/, $input));
	}
}
##############################

##############################
# gets the database (pattern) information of a file name
# input: mapped or unmapped file
# output: "_".pattern
sub get_pattern_split_files {
	my $file = $_[0];
	my $file_extension = $_[1];
	
	my $basename = fileparse($file, $file_extension);
	my @split = split(/_/, $basename);
	return "_".$split[-1];
}

##############################
# gets the reads that occur in one file but not in a array of files
# input: path of output directory, file with all reads, array of files to subtract reads, output file extension
# output: error message or 1
sub get_survived_reads {
	#################
	# save input parameters
	my $output_dir = $_[0];
	my $file = $_[1];
	my @files_removed = @{$_[2]};
	my $ofile_extension = $_[3];
	#################
	
	my @suffix = (".rc", ".map", "fasta", ".fa");
	
	# create output directory
	my $res = create_directory($output_dir); return $res if ($res ne 1);
	
	my $basename = fileparse($file, @suffix);
	
	my %reads = ();
	
	#################
	# save read and sequence information about combined file
	# open file handler
	open (IN, "<", $file) or return "Cannot open file $file!\n";
	while (<IN>) {
		my $line = $_;
		$line =~ tr/[\r\n]//d; # remove \n and \r characters
		if ($line =~ /^>(.*) count=\d+/) {
			my $id = $1;
			my $seq = <IN>;
			$seq =~ tr/[\r\n]//d; # remove \n and \r characters
			$reads{$id} = [$line, $seq];
		}
	}	
	# close file handler
	close IN;
	#################
	
	#################
	# delete reads occuring in files
	foreach my $remove_file (@files_removed) {
		open (IN, "<", $remove_file) or return "Cannot open file $remove_file!\n";
		while (<IN>) {
			my $line = $_;
			$line =~ tr/[\r\n]//d; # remove \n and \r characters
			if ($line =~ /^>(.*) count=\d+/) {
				my $id = $1;
				if (exists $reads{$id}) {
					delete $reads{$id};
				}
			}
		}
		close IN;
	}
	#################
	
	#################
	# print survived reads
	open (OUT, ">", $output_dir.$basename.$ofile_extension) or return "Cannot create file $output_dir$basename$ofile_extension!\n";
	foreach my $read (keys %reads) {
		print OUT $reads{$read}[0]."\n";
		print OUT $reads{$read}[1]."\n";
	}
	close OUT;
	#################
	
	return 1;
}
##############################

##############################
# prints used parameters to parameters.txt in pipeline output directory
# input: NA
# output: error message or 1
sub print_used_parameters {
	open (IN, "<", $config_file) or return "Cannot open file $config_file!\n";
	open (OUT, ">", $output_dir."parameters.txt") or return "Cannot create file ${output_dir}parameters.txt!\n";
	print OUT "######################################################\n";
	print OUT "######## PARAMETER SETTINGS FOR RUN $current_time ########\n";
	print OUT "######################################################\n";
	print OUT "##### Specified Input #####\n";
	($input ne "") ? print OUT "input: ".$input."\n" : print OUT "input: - \n";
	($output_dir ne "") ? print OUT "dir: ".$output_dir."\n" : print OUT "dir: - \n";
	print OUT "files: @input_files\n";
	my @settings = <IN>;
	for (my $i = 3; $i < scalar(@settings); $i++) {
		print OUT $settings[$i];
	}
	close IN;
	close OUT;
	
	return 1;
}
##############################

##############################
# quits the programm with a error message
# input: NA
# output: NA
sub quit {
	print "--------------------------\n";
	print "----- error detected -----\n";
	print "--------------------------\n";
	exit;
}
##############################

##############################
# saves map/unmap and sequence information in a hash
# input: mapped or unmapped file
# output: hash of mapped information and sequence
sub read_file_split_files {
	my $file = $_[0];
	my %hash_seq = ();
	open (IN, "<", $file) or check_status("Cannot open file $file!\n");
	while(<IN>) {
		my $line = $_;
		$line =~ tr/[\r\n]//d; # remove \n and \r characters
		if ($line =~ /^>read\d+\scount=\d+(.*)/) {
			my $tmp = $1;
			my $seq = <IN>;
			$seq =~ tr/[\r\n]//d; # remove \n and \r characters
			$hash_seq{$seq} = $tmp;
		}
	}
	close IN;
	
	return \%hash_seq;
}
##############################

##############################
# splits the reads of the combined mapped and unmapped files to each sample file according to mapping table
# input: mapping table file, array of files mapped to filter databases, array of files unmapped to filter databases
# array of files mapped to known databases, array of files unmapped to known databases, survived known reads, survived filter reads
# output: NA
sub split_files {
	#################
	# save input parameters
	my $read_mapping_file = $_[0];
	my @filter_files_map = @{$_[1]};
	my @filter_files_unmap = @{$_[2]};
	my @known_files_map = @{$_[3]};
	my @known_files_unmap = @{$_[4]};
	my $survived_known_file = $_[5];
	my $survived_filter_file = $_[6];
	#################
	
	my $dir_filtered = $output_dir."filtered/";
	my $dir_known = $output_dir."known/";
	
	#################
	# create directories
	if (scalar(@filter_files_map) > 0) {
		create_directory($dir_filtered);
		create_directory($dir_filtered."mapped/");
		create_directory($dir_filtered."unmapped/");
		create_directory($dir_filtered."survived/");
	}
	if (scalar(@known_files_map) > 0) {
		create_directory($dir_known);
		create_directory($dir_known."mapped/");
		create_directory($dir_known."unmapped/");
		create_directory($dir_known."survived/");
	}
	#################
	
	my @basenames = ();
	
	#################
	# save information of combined read mapping file
	my %hash_combine = ();
	open (IN, "<", $read_mapping_file) or check_status("Cannot open file $read_mapping_file!\n");
	while (<IN>) {
		my $line = $_;
		$line =~ tr/[\r\n]//d; # remove \n and \r characters
		my @split = split(/\t/, $line);
		if ($line !~ /^[ACGTUNacgtun]/) {
			foreach my $s (@split) {
				push (@basenames, $s) if ($s ne "");
			}
		} else {
			my $seq = $split[0];
			my $j = 0;
			for (my $i = 1; $i < scalar(@split)-1; $i=$i+2) {
				$hash_combine{$seq}{$basenames[$j]} = [$split[$i], $split[$i+1]];
				$j++;
			}
		}
	}
	close IN;
	#################
	
	#################
	# split and write each combined file to sample files
	foreach my $file (@filter_files_map) {		
		my $hash_seq = read_file_split_files($file);
		write_split_files($dir_filtered."mapped/", \@basenames, get_pattern_split_files($file, ".map"), \%hash_combine, $hash_seq, ".map");
	}
	foreach my $file (@filter_files_unmap) {		
		my $hash_seq = read_file_split_files($file);
		write_split_files($dir_filtered."unmapped/", \@basenames, get_pattern_split_files($file, ".rc"), \%hash_combine, $hash_seq, ".rc");
	}
	foreach my $file (@known_files_map) {		
		my $hash_seq = read_file_split_files($file);
		write_split_files($dir_known."mapped/", \@basenames, get_pattern_split_files($file, ".map"), \%hash_combine, $hash_seq, ".map");
	}
	foreach my $file (@known_files_unmap) {
		my $hash_seq = read_file_split_files($file);
		write_split_files($dir_known."unmapped/", \@basenames, get_pattern_split_files($file, ".rc"), \%hash_combine, $hash_seq, ".rc");
	}
	if ($survived_known_file ne "") {
		my $hash_seq = read_file_split_files($survived_known_file);
		write_split_files($dir_known."survived/", \@basenames, "", \%hash_combine, $hash_seq, ".map");
	}
	if ($survived_filter_file ne "") {
		my $hash_seq = read_file_split_files($survived_filter_file);
		write_split_files($dir_filtered."survived/", \@basenames, "_genome", \%hash_combine, $hash_seq, ".map");
	}
	#################
}
##############################

##############################
# checks all needed database files
# input: NA
# output: error message or 1
sub test_database_files {
	if (-d $path_species) {
		my $res = check_database_bowtie("genome", $path_species); return $res if ($res ne 1);
		
		if ($config->getProperty("filter_reads") eq "yes") {
			my $path_exclcriteria = $path_species."exclusioncriteria/";
			if (-d $path_exclcriteria) {
				# get all possible filter criteria
				opendir (DIR, $path_exclcriteria) or return "Cannot open directory $path_exclcriteria\n";
				# get directory names of all filter criteria
				my @filter_dirs = readdir(DIR);
				# loop over each filter criteria
				foreach my $dir_filter_name (@filter_dirs) {
					if ($config->getProperty($dir_filter_name, "") eq "yes") {
						$res = check_database_bowtie($dir_filter_name, $path_exclcriteria); return $res if ($res ne 1);
					}
				}
			} else {
				return "Cannot find exclusioncriteria folder!\n";
			}
		}
		
		if ($config->getProperty("detect_known") eq "yes") {
			my $path_microRNA = $path_species."microRNA/";
			if (-d $path_microRNA) {
				$res = check_database_bowtie("mature", $path_microRNA); return $res if ($res ne 1);
				$res = check_database_bowtie("precursors", $path_microRNA); return $res if ($res ne 1);
			} else {
				return "Cannot find microRNA folder!\n";
			}
		}
		
		if ($config->getProperty("predict_target") eq "yes") {
			$res = check_database("mRNA", $path_species); return $res if ($res ne 1);
		}
		
		return 1;
	} else {
		return "Cannot find data folder for species $species!\n";
	}
}
##############################	

##############################
# tests and returns the format of given input files
# input: NA
# output: error message or "rc" or "fastq"
sub test_input_files {
	my $format = "";
	foreach my $f (@input_files) {
		print "      \t".$f.": ";
		my $file = fileparse($f);	
		open (IN, "<", $f) or return "Cannot open file $file!\n";
		my $cnt = 0;
		while (my $header = <IN>) {
			$cnt++;
			if ($header =~ /^@/) {
				if ($format eq "rc" && $cnt == 1) {
					return "Use only one file format (FASTQ or read-count) for all files!\n";
				} elsif ($format eq "rc" && $cnt > 1) {
					return "Line $cnt in file $file should start with '>'!\n";
				}
				$format = "fastq";
				my $sequence = <IN>;
				$sequence =~ tr/[\r\n]//d; # remove \n and \r characters
				$cnt++;
				if ($sequence !~ /^[ACGTNacgtn]+$/) {
					return "The sequence in line $cnt in file $file contains characters other than A, C, G, T and N!\n";
				} else {
					my $header2 = <IN>;
					$cnt++;
					if ($header2 !~ /^\+/) {
						return "Line $cnt in file $file should start with '\+' \n";
					}
					my $quality = <IN>;
					$cnt++;
				}	
			} elsif ($header =~ /^>/) {
				if ($format eq "fastq" && $cnt == 1) {
					return "Use only one file format (FASTQ or read-count) for all files!\n";
				} elsif ($format eq "fastq" && $cnt > 1) {
					return "Line $cnt in file $file should start with '\@'!\n";
				}
				$format = "rc";
				
				if ($header !~ />.*count=\d+/) {
					return "Line $cnt in file $file doesn't match a header of a read-count file!\n";
				}
				my $sequence = <IN>;
				$cnt++;
				$sequence =~ tr/[\r\n]//d; # remove \n and \r characters
				if ($sequence !~ /^[ACGTNacgtn]+$/) {
					return "The sequence in line $cnt in file $file contains characters other than A, C, G, T and N!\n";
				}
			} else {
				return "File $file is not a FASTQ or read-count file at line $cnt!\n";
			}
		}
		close IN;
		print "done\n";
	}
	return $format;
}		
##############################

##############################
# writes map files per sample
# input: path of output directory, array of basenames, pattern, hash of combined reads, hash of sequences, output file extension
# output: NA
sub write_split_files {
	#################
	# save input parameters
	my $output_dir = $_[0];
	my @basenames = @{$_[1]};
	my $pattern = $_[2];
	my %hash_combine = %{$_[3]};
	my %hash_seq = %{$_[4]};
	my $file_extension = $_[5];
	#################

	foreach my $sample (@basenames) {	
		open (OUTMAP, ">", $output_dir.$sample.$pattern.$file_extension) or check_status("Cannot create file $output_dir$sample$pattern$file_extension!\n");
		foreach my $seq (keys %hash_seq) {
			if ($hash_combine{$seq}{$sample}[0] ne "-") {
				print OUTMAP ">".$hash_combine{$seq}{$sample}[0];
				print OUTMAP " count=".$hash_combine{$seq}{$sample}[1];
				print OUTMAP $hash_seq{$seq}."\n";
				print OUTMAP $seq."\n";
			}
		}
		close OUTMAP;
	}
}
##############################
