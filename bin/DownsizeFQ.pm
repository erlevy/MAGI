#!/usr/bin/perl
#03/29/2012
#Author: Petra Stepanowsky

package DownsizeFQ;

use warnings;
use strict;
use File::Basename;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);


use Time::HiRes qw(gettimeofday tv_interval);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = qw(downsize);
@EXPORT_OK   = qw(downsize);
%EXPORT_TAGS = (DEFAULT => [qw(&downsize)]);

my $num_Ns_total = "";
my $offset = "";
my $qual_type = "";
my $qual_threshold = "";
my $adapter3 = "";
my $adapter_proportion = "";
my $num_mismatches = "";
my $flag_count_mis_N = "";
my $min_length = "";
my $max_length = "";
my $summary = "";
my $qualfile = "";
my $trim = "";
my @adapter3_chars = ();

# maximum quality value to achieve
my $max_score = 60;

##############################
# creates a directory
# input: directory name/path
# output: error message or 1
sub create_directory {
	unless(-d $_[0]){
		mkdir $_[0] or return "Cannot create directory ($_[0])\n";
	}
	return 1;
}
##############################

##############################
# collapses, filters and trims 3' adapter reads
# input: array of files in FASTQ-format, path of output directory, number of Ns in a sequence, 
# offset of FASTQ format, type of quality function (min, mean, median), threshold of quality,
# 3' adapter sequence, minimum number of nts of 3' adapter that have to match with read sequence,
# maximum number of allowed mismatches, flag to count mismatch with N or not, minimum length,
# maximum length, flag to create summary file, flag to create quality file, flag to trim reads
# output: error message or 1
sub downsize {	
	#################
	# save input parameters
	my @input_files = @{$_[0]};
	my $output_dir = $_[1]."/downsized/";
	$num_Ns_total = $_[2];
	$offset = $_[3];
	$qual_type = $_[4];
	$qual_threshold = $_[5];
	$adapter3 = $_[6];
	$adapter_proportion = $_[7];
	$num_mismatches = $_[8];
	$flag_count_mis_N = $_[9];
	$min_length = $_[10];
	$max_length = $_[11];
	$summary = $_[12];
	$qualfile = $_[13];
	$trim = $_[14];
	#################
	
	my $qual_dir = $output_dir."qual/";
	my $rc_dir = $output_dir."rc/";
	
	# split 3' adapter sequence
	@adapter3_chars = split(//, $adapter3); # characters of the 3' adapter sequence

	my $TOTAL = "all";
	my $PREPROCESSED = "more_than_".($num_Ns_total)."_N";
	my $HIGHQUAL = "high_quality";
	my $LOWQUAL = "low_quality";
	my $NOADAPTER = "no_adapter";
	my $REMAINDER = "remainder";
	my $TOOSHORTREADS = "shorter_".$min_length;
	my $TOOLONGREADS = "longer_".$max_length;
	
	my $summary_total_file = $output_dir."downsize_summary_total.txt";
	my $summary_unique_file = $output_dir."downsize_summary_unique.txt";

	my $res = create_directory($output_dir); return $res if ($res ne 1);
	$res = create_directory($qual_dir); return $res if ($res ne 1);
	$res = create_directory($rc_dir); return $res if ($res ne 1);

	# number of all files
	my $num_files = scalar(@input_files);
	# count processed files
	my $cnt_files = 1;
	# file extensions for FASTQ files
	my @suffix = (".fastq", ".fq");

	(open (SUMTOTAL, ">", $summary_total_file) or return "Cannot create file $summary_total_file!\n") if ($summary eq "yes");
	(open (SUMUNIQUE, ">", $summary_unique_file) or return "Cannot create file $summary_unique_file!\n") if ($summary eq "yes");

	if ($summary eq "yes") {
		print SUMTOTAL "sample\t$TOTAL\t$PREPROCESSED\t$HIGHQUAL\t$LOWQUAL\t$NOADAPTER\t$TOOSHORTREADS\t$TOOLONGREADS\t$REMAINDER\n";
		print SUMUNIQUE "sample\t$TOTAL\t$PREPROCESSED\t$HIGHQUAL\t$LOWQUAL\t$NOADAPTER\t$TOOSHORTREADS\t$TOOLONGREADS\t$REMAINDER\n";
	}

	foreach my $file (@input_files) {
	
		# file name without extension
		my $filebasename = fileparse($file, @suffix);
	
		#################
		# print current processing file
		print ">".$filebasename." (".$cnt_files."/".$num_files.")";
		print "\n\tcollapsing reads";
		#################
		
		# open filehandler
		open (IN, "<", $file) or return "Cannot open file $file!\n";
		
		# hash to save unique reads
		# key: read sequence
		# value: array containing two elements, index 0: number of high quality lines
		#                                       index 1: number of low quality lines
		my %reads = ();
		
		#################
		# number of all reads in a file
		my $all_total = 0;
		my $all_unique = 0;
		#################
		
		#################
		# number of reads containing too many 'N' characters
		my $discardedN_unique = 0;
		my $discardedN_total = 0;
		#################
		
		#################
		# number of 'low quality' reads
		my $lowqual_unique = 0;
		my $lowqual_total = 0;
		#################
		
		#################
		# number of 'high quality' reads
		my $highqual_unique = 0;
		my $highqual_total = 0;
		#################
		
		#################
		# number of 'high quality' reads with no found adapter
		my $noadapter_unique = 0;
		my $noadapter_total = 0;
		#################
		
		#################
		# number of trimmed 'high quality' reads that are too short
		my $short_unique = 0;
		my $short_total  = 0;
		#################
		
		#################
		# number of trimmed 'high quality' reads that are too long
		my $long_unique = 0;
		my $long_total  = 0;
		#################
		
		#################
		# number of trimmed 'high quality' reads that are written to an output file
		my $remainder_unique = 0;
		my $remainder_total = 0;
		#################
		
		#################
		# count variables for quality file
		my @min_quality = ();
		my @max_quality = ();
		my @mean_quality = ();
		my @median_qual_all = ();
		my @median_quality = ();
		my @thirdQ = ();
		my @firstQ = ();
		my @sd_quality = ();
		
		# array index 0: total count, index 1: unique count
		my @num_A = ();
		my @num_C = ();
		my @num_G = ();
		my @num_T = ();
		my @num_N = ();
		#################
		
		#################
		# count sequence length of 'high quality' reads after 3' adapter trimming
		my @read_lengths = ();
		#################
		
		#################
		# read in file line by line
		while (my $header = <IN>) { # '@' line
		
			if ($header =~ /^@/) { # first line must start with '@'
				
				#################
				# read sequence
				my $sequence = <IN>;
				$sequence =~ tr/[\r\n]//d; # remove \n and \r characters
				$sequence = uc($sequence);
				#################
				
				if ($sequence =~ /^[ACGTNacgtn]+$/) { # sequence only contains A, C, G, T or N

					# initialize array in hash if sequence doesn't exist
					$reads{$sequence} = [0,0] unless (exists $reads{$sequence});
					
					#################
					# '+' line
					my $header2 = <IN>;
					#################
					
					if ($header2 =~ /^\+/) {
					
						#################
						# quality line
						my $quality = <IN>;
						$quality =~ tr/[\r\n]//d; # remove \n and \r characters
						#################
						
						# values of quality line
						my @quality_values = ();
						
						#################
						# calculate quality distribution per position
						# split quality line
						my $qual_ind = length($quality) - 1;

						while (my $q = ord(chop $quality)) {
							# convert quality character to ascii and subtract the offset
							my $qual_score = $q - $offset;
							
							#################
							# calculate values for quality file
							if ($qualfile eq "yes") { # if quality file is desired
								#################
								# sum up quality values to calculate mean quality
								$mean_quality[$qual_ind] += $qual_score;
								#################
							
								#################
								# initialize first minimum and maximum value
								if($all_total eq 0) {
									$min_quality[$qual_ind] = 1000;
									$max_quality[$qual_ind] = -1000;
									@{$median_qual_all[$qual_ind]} = (0) x $max_score;
								}
								#################
								
								#################
								# set minimum quality
								if ($min_quality[$qual_ind] > $qual_score) {
									$min_quality[$qual_ind] = $qual_score;
								}
								#################
								
								#################
								# set maximum quality
								if ($max_quality[$qual_ind] < $qual_score) {
									$max_quality[$qual_ind] = $qual_score;
								}
								#################
								
								#################
								# set median quality
								$median_qual_all[$qual_ind][$qual_score]++;
								#################
							}
							#################
							
							$quality_values[$qual_ind] = $qual_score;
							--$qual_ind;
						}
						#################
						
						#################
						# check if read is a 'high quality' read or not
						if (has_high_quality(\@quality_values)) {
							$reads{$sequence}[0]++;	
						} else {
							$reads{$sequence}[1]++;
						}
						#################
						
						# count reads in the file
						$all_total++;
					}	
				}
			}
		}
		#################
		
		# close filehandler
		close IN;
		
		# number of all unique reads
		$all_unique = scalar keys %reads;
		
		#################
		# calculate mean quality per position
		if ($qualfile eq "yes") { # if quality file is desired
			foreach my $mean (@mean_quality) {
				# calculate mean for each position
				$mean = $mean / $all_total;
			}
			
			#################
			# for each position calculate standard deviation, median
			for (my $k = 0; $k < scalar(@median_qual_all); $k++) {
				my @tmp = @{$median_qual_all[$k]};
				my $tmp_size = scalar(@tmp);

				my $s = 0;
				$s += $_ for @tmp;
				
				#################
				#  calculate standard deviation for each position
				my $sum = 0;
				for (my $l = 0; $l < $tmp_size; $l++) {
					$sum += (($l - $mean_quality[$k]) ** 2) * $tmp[$l];
				}
				$sd_quality[$k] = sqrt($sum / $s);
				#################
					
				#################
				# calculate median quality, upper and lower quartile for each position
				$median_quality[$k] = calc_median(0, $tmp_size - 1, \@tmp);
				$firstQ[$k] = calc_firstQ(0, $tmp_size - 1, \@tmp);
				$thirdQ[$k] = calc_thirdQ(0, $tmp_size - 1, \@tmp);
				#################
			}
			#################
		}
		#################
		
		# hash to save unique trimmed reads
		# key: trimmed read sequence
		# value: read count
		my %remainder_reads = ();
		
		print "\n\tprocessing unique reads (".$all_unique." reads)";
					
		my $first_flag = 1;
		#################
		# process set of unique reads
		foreach my $read_seq (keys %reads) {
			
			# expression level (with high and low quality lines)
			my $read_count = $reads{$read_seq}[0] + $reads{$read_seq}[1];

			#################
			# calculate base distribution per position
			if ($qualfile eq "yes") { # if quality file is desired
				my $tmp = $read_seq;
				my $ind = length($tmp) - 1 ;
				
				while (my $c = chop $tmp) {
					#################
					# initialize first A, C, G, T and N
					if($first_flag) {
						$num_A[$ind][0] = 0;
						$num_A[$ind][1] = 0;
						$num_C[$ind][0] = 0;
						$num_C[$ind][1] = 0;
						$num_G[$ind][0] = 0;
						$num_G[$ind][1] = 0;
						$num_T[$ind][0] = 0;
						$num_T[$ind][1] = 0;
						$num_N[$ind][0] = 0;
						$num_N[$ind][1] = 0;
					}
					#################

					if ($c eq "A") {
						$num_A[$ind][0] += $read_count;
						$num_A[$ind][1]++;
					} elsif ($c eq "C") {
						$num_C[$ind][0] += $read_count;
						$num_C[$ind][1]++;
					} elsif ($c eq "G") {
						$num_G[$ind][0] += $read_count;
						$num_G[$ind][1]++;
					} elsif ($c eq "T") {
						$num_T[$ind][0] += $read_count;
						$num_T[$ind][1]++;
					} elsif ($c eq "N") {
						$num_N[$ind][0] += $read_count;
						$num_N[$ind][1]++;
					}	
					
					--$ind;
				}
			}
			#################
			
			#################
			# count character 'N' in read sequence
			if ($read_seq =~ tr/N// > $num_Ns_total) {
				# number of N's in read sequence > number of N's specified
				$discardedN_unique++;
				$discardedN_total += $read_count;
			} else {
				#################
				# sum up low quality lines of current read
				if ($reads{$read_seq}[1] > 0) {
					# number of low quality lines greater than 0
					$lowqual_unique++;
					$lowqual_total += $reads{$read_seq}[1];
				}
				#################
			
				#################
				# process 'high quality' read
				if ($reads{$read_seq}[0] > 0) {
					# number of high quality lines greater than 0
						
					#################
					# sum up high quality lines of current read
					$highqual_unique++;
					$highqual_total += $reads{$read_seq}[0];
					#################
						
					if ($trim eq "yes") {
						#################
						# trim 3' adapter
						my $trimmed_read_seq = remove_adapter($read_seq);
						#################
						
						#################
						# count length of trimmed read
						$read_lengths[length($trimmed_read_seq)][0]++;
						$read_lengths[length($trimmed_read_seq)][1] += $reads{$read_seq}[0];
						#################
						
						#################
						# no 3' adapter sequence was found
						if ($trimmed_read_seq eq $read_seq) {
							$noadapter_unique++;
							$noadapter_total += $reads{$read_seq}[0];
						}
						#################
						
						#################
						# filter too long and too short reads
						if (length($trimmed_read_seq) < $min_length) {
							# trimmed read sequence length is too short
							$short_unique++;
							$short_total += $reads{$read_seq}[0];
						} elsif (length($trimmed_read_seq) > $max_length) {
							# trimmed read sequence length is too long
							$long_unique++;
							$long_total += $reads{$read_seq}[0];
						} else {
							# save trimmed read
							$remainder_reads{$trimmed_read_seq} += $reads{$read_seq}[0];
						}
						#################
					} else {
						$remainder_reads{$read_seq} += $reads{$read_seq}[0];

						#################
						# count length of trimmed read
						$read_lengths[length($read_seq)][0]++;
						$read_lengths[length($read_seq)][1] += $reads{$read_seq}[0];
						#################
					}
				}
			#################
			}
			#################
			
			$first_flag = 0;
		}
		#################
		
		#################
		# write trimmed remainder reads to output file
		# open file handler
		open (OUT, ">", $rc_dir.$filebasename.".rc") or return "Cannot create file $rc_dir$filebasename.rc!\n";

		print "\n\twriting remainder reads (".(scalar keys %remainder_reads)." reads)";

		# count trimmed reads, add number of processed read to generate unique identifier
		my $cnt_reads = 1;

		foreach my $read_seq (sort {$remainder_reads{$b} <=> $remainder_reads{$a}} keys %remainder_reads) {
			$remainder_unique++;
			$remainder_total += $remainder_reads{$read_seq};

			print OUT ">read".$cnt_reads." count=".$remainder_reads{$read_seq}."\n";
			print OUT $read_seq."\n";

			$cnt_reads++;
		}
		# close file handler
		close OUT;
		#################
			
		if ($summary eq "yes") {
			print "\n\twriting summary";

			print SUMTOTAL $filebasename."\t".$all_total."\t";
			print SUMUNIQUE $filebasename."\t".$all_unique."\t";
			print SUMTOTAL $discardedN_total."\t";
			print SUMUNIQUE $discardedN_unique."\t";
			print SUMTOTAL $highqual_total."\t";
			print SUMUNIQUE $highqual_unique."\t";
			print SUMTOTAL $lowqual_total."\t";
			print SUMUNIQUE $lowqual_unique."\t";
			($trim eq "yes") ? print SUMTOTAL $noadapter_total."\t" : print SUM "NA\tNA\t";
			($trim eq "yes") ? print SUMUNIQUE $noadapter_unique."\t" : print SUM "NA\tNA\t";
			($trim eq "yes") ? print SUMTOTAL $short_total."\t" : print SUM "NA\tNA\t";
			($trim eq "yes") ? print SUMUNIQUE $short_unique."\t" : print SUM "NA\tNA\t";
			($trim eq "yes") ? print SUMTOTAL $long_total."\t" : print SUM "NA\tNA\t";
			($trim eq "yes") ? print SUMUNIQUE $long_unique."\t" : print SUM "NA\tNA\t";
			print SUMTOTAL $remainder_total;
			print SUMUNIQUE $remainder_unique;

			print SUMTOTAL "\n" if ($cnt_files < $num_files);
			print SUMUNIQUE "\n" if ($cnt_files < $num_files);
		}

		if ($qualfile eq "yes") {
			print "\n\twriting quality file";

			my $res = print_base_quality($qual_dir.$filebasename."_basequality.txt", \@mean_quality, \@sd_quality, \@min_quality, \@thirdQ, \@median_quality, \@firstQ, \@max_quality,
								\@num_A, \@num_C, \@num_G, \@num_T, \@num_N, \@median_qual_all);
			return $res if ($res ne 1);
			$res = print_read_lengths($qual_dir.$filebasename."_readlength.txt", \@read_lengths);
			return $res if ($res ne 1);
		}
		print "\n";
		$cnt_files++;
	}
	# close file handler
	close SUMTOTAL if ($summary eq "yes");
	close SUMUNIQUE if ($summary eq "yes");

	return 1;
}
##############################

	
##############################
# checks if a given array containing quality values, contains 'high quality' values
# input: array containing quality values
# output: if the quality is higher than a certain threshold, returns 1 (true), otherwise 0 (false)
sub has_high_quality {
	if ($qual_type eq "min") {
		return &min_quality($_[0]);
	} elsif ($qual_type eq "median") {
		return &median_quality($_[0]);
	} elsif ($qual_type eq "mean") {
		return &mean_quality($_[0]);
	}
}
##############################
	
##############################
# checks the mean quality of a given array containing quality values
# input: array containing quality values
# output: if the quality is higher than a certain threshold, returns 1 (true), else 0 (false)
sub mean_quality {
	#################
	# save input parameter
	my @values = @{$_[0]};
	#################
	
	my $sum = 0;

	#################
	# calculate sum of given values
	foreach my $value (@values){
		$sum += $value;
	}
	#################

	# calculate mean of given values
	my $mean = $sum / scalar(@values);

	#################
	# if mean is smaller than a certain threshold, return 1, else return 0
	if ($mean < $qual_threshold) {
		# mean is smaller than threshold
		return 0;
	} else {
		# mean is equal or greater than threshold
		return 1;
	}
	#################
}
##############################

##############################
# calculates the lower quartile
# input: index left, index right, array
# output: lower quartile
sub calc_thirdQ {
	my $i = $_[0];
	my $j = $_[1];
	my @y = @{$_[2]};

	while ($j > $i && !(($y[$i] == 0 && $y[$j] == 0) && $i + 1 >= $j)) {
		if ($y[$i] > 0) {
			--$y[$i];
		} else {
			while ($y[$i] == 0 && $i < $j) {
				$i++;
			}
			if ($i < $j) {
				--$y[$i];
			}
		}

		for (my $k = 0; $k < 3; $k++) {
			if ($y[$j] > 0) {
				--$y[$j];
			} else {
				while ($y[$j] == 0 && $j > $i) {
					$j--;
				}
				if ($j > $i) {
					--$y[$j];
				}
			}
		}
	}
	($i == $j) ? return $i : (($i + $j) / 2);
}
##############################

##############################
# calculates the median
# input: index left, index right, array
# output: median
sub calc_median {
	my $i = $_[0];
	my $j = $_[1];
	my @y = @{$_[2]};

	while ($j > $i && !(($y[$i] == 0 && $y[$j] == 0) && $i + 1 >= $j)) {
		if ($y[$i] > 0) {
			--$y[$i];
		} else {
			while ($y[$i] == 0 && $i < $j) {
				$i++;
			}
			if ($i < $j) {
				--$y[$i];
			}
		}
		
		if ($y[$j] > 0) {
			--$y[$j];
		} else {
			while ($y[$j] == 0 && $j > $i) {
				$j--;
			}
			if ($j > $i) {
				--$y[$j];
			}
		}
	}
	($i == $j) ? return $i : (($i + $j) / 2);
}
##############################

##############################
# calculates the upper quartile
# input: index left, index right, array
# output: upper quartile
sub calc_firstQ {
	my $i = $_[0];
	my $j = $_[1];
	my @y = @{$_[2]};
	
	while ($j > $i && !(($y[$i] == 0 && $y[$j] == 0) && $i + 1 >= $j)) {
		for (my $k = 0; $k < 3; $k++) {
			if ($y[$i] > 0) {
				--$y[$i];
			} else {
				while ($y[$i] == 0 && $i < $j) {
					$i++;
				}
				if ($i < $j) {
					--$y[$i];
				}
			}
		}
		
		if ($y[$j] > 0) {
			--$y[$j];
		} else {
			while ($y[$j] == 0 && $j > $i) {
				$j--;
			}
			if ($j > $i) {
				--$y[$j];
			}
		}
	}
	($i == $j) ? return $i : (($i + $j) / 2);
}
##############################
	
##############################
# checks the median quality of a given array containing quality values
# input: array containing quality values
# output: if the quality is higher than a certain threshold, returns 1 (true), else 0 (false)
sub median_quality {
	#################
	# save input parameter
	my @values = @{$_[0]};
	#################
	
	my $median = 0;
	
	# sort values
	@values = sort {$a <=> $b} @values;

	# index of middle element
	my $middle = scalar(@values) / 2;
	#################
	# if the number of elements is odd, the median is the middle element,
	# otherwise the median is the mean of the two middle elements
	if (scalar(@values) % 2) {
		# number of elements is odd
		$median = $values[$middle];
	} else {
		# number of elements is even
		# median is the mean of the two middle elements
		$median = ($values[$middle] + $values[$middle - 1]) / 2;
	}
	#################
	
	#################
	# if median is smaller than a certain threshold, return 1, else return 0
	if ($median < $qual_threshold) {
		# median is smaller than threshold
		return 0;
	} else {
		# median is equal or greater than threshold
		return 1;
	}
	#################
}
##############################

	
##############################
# checks the minimum quality of a given array containing quality values
# input: array containing quality values
# output: if the quality is higher than a certain threshold, returns 1 (true), else 0 (false)
sub min_quality {
	#################
	# save input parameter
	my @values = @{$_[0]};
	#################
	
	# sort values
	@values = sort {$a <=> $b} @values;

	#################
	# if the first element of the sorted array is smaller than a threshold,
	# the minimum of the values is smaller than a threshold, return 1, otherwise return 0
	if ($values[0] < $qual_threshold) {
		# first element is smaller than threshold
		return 0;
	} else {
		# first element is equal or greater than threshold
		return 1;
	}
	#################
}
##############################	


##############################
# writes the base and quality distribution to a file
# input: name of output file, array containing mean quality per position, minimum quality per
# position, maximum quality per position, number of A's, number of C's, number of G's, number
# of T's and number of N's, respectively
# output: error message or 1
sub print_base_quality {
	#################
	# save input parameters
	my $file = $_[0];
	my @mean = @{$_[1]};
	my @sd = @{$_[2]};
	my @min = @{$_[3]};
	my @firstQ = @{$_[4]};
	my @median = @{$_[5]};
	my @thirdQ = @{$_[6]};
	my @max = @{$_[7]};
	my @num_A = @{$_[8]};
	my @num_C = @{$_[9]};
	my @num_G = @{$_[10]};
	my @num_T = @{$_[11]};
	my @num_N = @{$_[12]};
	my @qual_all = @{$_[13]};
	#################

	# open file handler	
	open(OUT, ">", $file) or return "Cannot create file $file!\n";

	# print header
	print OUT "position\tmean\tsd\tmin\t3Q\tmedian\t1Q\tmax\tGCratio_total\tGCratio_unique\tA_total\tC_total\tG_total\tT_total\tN_total\tA_unique\tC_unique\tG_unique\tT_unique\tN_unique\tcount_per_quality\n";

	#################
	# foreach position print mean, min, max, GC ratio and number of each nucleotide
	for (my $i = 0; $i <= $#mean; $i++) {
		# position
		print OUT ($i+1)."\t";
		
		# mean quality per position
		print OUT $mean[$i]."\t";
		
		# standard deviation quality per position
		print OUT $sd[$i]."\t";
		
		# minimum quality per position
		print OUT $min[$i]."\t";
		
		# lower quartile per position
		print OUT $firstQ[$i]."\t";
		
		# median quality per position
		print OUT $median[$i]."\t";
		
		# upper quartile per position
		print OUT $thirdQ[$i]."\t";
		
		# maximum quality per position
		print OUT $max[$i]."\t";
		
		# GC ratio per position total
		print OUT (($num_G[$i][0] + $num_C[$i][0]) / ($num_A[$i][0] + $num_C[$i][0] + $num_G[$i][0] + $num_T[$i][0] + $num_N[$i][0]))."\t";
		# GC ratio per position unique
		print OUT (($num_G[$i][1] + $num_C[$i][1]) / ($num_A[$i][1] + $num_C[$i][1] + $num_G[$i][1] + $num_T[$i][1] + $num_N[$i][1]))."\t";
		
		#number of bases per position total
		print OUT $num_A[$i][0]."\t";
		print OUT $num_C[$i][0]."\t";
		print OUT $num_G[$i][0]."\t";
		print OUT $num_T[$i][0]."\t";
		print OUT $num_N[$i][0]."\t";
		
		# number of base per position unique
		print OUT $num_A[$i][1]."\t";
		print OUT $num_C[$i][1]."\t";
		print OUT $num_G[$i][1]."\t";
		print OUT $num_T[$i][1]."\t";
		print OUT $num_N[$i][1]."\t";
		
		print OUT "@{$qual_all[$i]}\n";
	}
	#################

	# close file handler
	close OUT;
	
	return 1;
}
##############################
	
##############################
# writes read lengths and count numbers to a file
# input: name of output file, array containing read lengths and count numbers
# output: error message or 1
sub print_read_lengths {
	#################
	# save input parameters
	my $file = $_[0];
	my @read_lengths = @{$_[1]};
	#################
	
	# open file handler
	open (OUT, ">", $file) or return "Cannot create file $file!\n";

	# print header
	print OUT "read_length\tunique\ttotal\n";

	#################
	# for each read length print length and total and unique count number
	for (my $i = 0; $i <= $#read_lengths; $i++) {
		if (defined $read_lengths[$i]) {
			# read length
			print OUT $i."\t";
			# number of unique 'high quality' reads
			print OUT $read_lengths[$i][0]."\t";
			# number of total 'high quality' reads
			print OUT $read_lengths[$i][1]."\n";
		} else {
			# no read with this certain read length exists
			print OUT $i."\t";
			print OUT "0\t0\n";
		}
	}
	#################
	
	# close file handler
	close OUT;
	
	return 1;
}
##############################

##############################
# removes the 3' adapter sequence of a given sequence
# input: sequence
# output: if a 3' adapter sequence is found, returns the trimmed sequence, otherwise returns the input sequence
sub remove_adapter {
	#################
	# save input parameter
	my $sequence = $_[0];
	#################
	
	# trimmed sequence
	my $sequence_trimmed = $sequence;

	#################
	# split input sequence
	# array containing characters
	my @sequence_chars = split(//, $sequence);
	#################

	#################
	# alignment of input sequence and 3' adapter sequence
	# outer loop for input sequence
	for (my $i = 0; $i < length($sequence); $i++) {

		#################
		# reset match and mismatch variable
		my $match = 0;
		my $mismatch = 0;
		#################

		#################
		# inner loop for 3' adapter sequence
		for (my $j = 0; $j < length($adapter3) && ($mismatch < $num_mismatches); $j++) {
			# last itertation if the end of the input sequence is reached
			last unless ($i + $j < length($sequence));

			#################
			# sum up match or mismatch
			if ($sequence_chars[$i + $j] eq $adapter3_chars[$j]) {
				# the characters of the two sequences are the same
				$match++;
			} else {
				# the characters of the two sequences are different
				if ($sequence_chars[$i + $j] ne "N" && $adapter3_chars[$j] ne "N") {
					# the characters are not a 'N'
					$mismatch++;
				} elsif (($sequence_chars[$i + $j] eq "N" || $adapter3_chars[$j] eq "N") && $flag_count_mis_N) {
					# at least one character is a 'N'
					# and the flag to track mismatches with a character 'N' is set
					$mismatch++;
				}
			}
			#################
		}
		#################

		if (($mismatch < $num_mismatches) && ($match >= $adapter_proportion)) {
			# the trimmed sequence is a substring of the input sequence
			# i is the index, where the 3' adapter sequence starts
			$sequence_trimmed = substr($sequence, 0, $i);

			last; # a possible alignment is found
		}
	}
	#################

	# return trimmed read
	return $sequence_trimmed;
}
##############################	

#loaded ok
1;