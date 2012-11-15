#!/usr/bin/perl
#created 03/29/2012, last updated 08/10/2012
#Author: Petra Stepanowsky

package Aligner;

use warnings;
use strict;
use File::Basename;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = qw(map);
@EXPORT_OK   = qw(map);
%EXPORT_TAGS = (DEFAULT  => [qw(&map)]);

my @suffix = (".fasta", ".fa", ".rc", ".sam", ".map");

##############################
# align given input file to given reference using Bowtie 2
# input: path to FASTA-formatted file, path to basename of Bowtie 2 index, path to output file, number of mismatches, Bowtie 2 parameter
# output: system call error code or 1
sub align {
	my $res = system("bowtie2 --gbar 32 -f $_[4] -N $_[3] -x $_[1] $_[0] -S $_[2] 2>/dev/null");
	$res == 0 ? return 1 : return $res;
}
##############################

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
# aligns given input files to a given reference using Bowtie 2 and splits the output files in mapped and unmapped files
# input: input files, path of basename of Bowtie 2 index, database name pattern, path of output directory, number of mismatches, optional: additional parameter for Bowtie 2 call
# output: error message or 1
sub map {
	#################
	# save input parameters
	# input files
	my @files = @{$_[0]};
	# basename of Bowtie 2 index
	my $reference = $_[1];
	# pattern of database; added to output file name
	my $pattern = $_[2];
	# output directory
	my $output_dir = $_[3];
	# number of mismatches
	my $mismatches = $_[4];
	# additional parameter for Bowtie 2
	my $parameter = "";
	if (exists $_[5]) { 
		$parameter = $_[5];
	}
	#################
	
	# array of output file names
	my @ofiles = ();
	
	#################
	# create alignment output directories
	my $res = create_directory($output_dir); return $res if ($res ne 1);
	$res = create_directory($output_dir."alignment/"); return $res if ($res ne 1);
	#################
	
	#################
	# handle file by file
	foreach my $f (@files) {
		my $basename = fileparse($f, @suffix);
		if ($basename =~ /\_/) {
			my @split = split(/\_/, $basename);
			if (scalar(@split) > 1) {
				$basename = $split[0]."_".$split[1];
			}
		}
		
		# remove first whitespace of files with 'count'-property in header
		my $f_tmp = remove_whitespaces_header($f, $output_dir); return $f_tmp if ($f_tmp !~ /.tmp$/);
		
		# output file name
		my $ofile = $output_dir."alignment/".$basename."_".$pattern.".sam";
		
		# align file to reference
		$res = align($f_tmp, $reference, $ofile, $mismatches, $parameter); return $res if ($res ne 1);
		
		# save output file name in array for further processing
		push(@ofiles, $ofile);
		# remove temporary generated file
		remove_file($f_tmp);
	}
	#################
	
	#################
	# create mapped and umapped directories
	$res = create_directory($output_dir."mapped/"); return $res if ($res ne 1);
	$res = create_directory($output_dir."unmapped/"); return $res if ($res ne 1);
	#################
	
	# split SAM-formatted files
	$res = split_SAM(\@ofiles, $output_dir."mapped/", $output_dir."unmapped/", $pattern); return $res if ($res ne 1);
	
	return 1;
}
##############################

##############################
# removes a file
# input: path to a file
sub remove_file {
	unlink $_[0];
}
##############################

##############################
# removes the first whitespace of the header in a read-count file
# input: path of read-count-formatted file, path of output directory
# output: error message or path of FASTA-formatted file with extension '.tmp'
sub remove_whitespaces_header {
	#################
	# save input parameters
	my $file = $_[0];
	my $output_dir = $_[1];
	#################
	
	#################
	# initialize variables
	my $basename = fileparse($file, @suffix);
	my $tmp_file = $output_dir."/".$basename.".tmp";
	#################
	
	# open input file handler
	open (IN, "<", $file) or return "Cannot open file $file!\n";
	# open output file handler
	open (OUT, ">", $tmp_file) or return "Cannot create file $tmp_file!\n";
	
	#################
	# read input file
	while (<IN>) {
		my $line = $_;
		$line =~ tr/[\r\n]//d; # remove \n and \r characters
		
		if ($line =~ /^>/) {
			if ($line =~ /^(>.*\scount=\d+)/) {
				$line = $1;
				$line =~ s/\s/\_/;
			} 
			print OUT $line."\n";
		} elsif ($line =~ /^[ACGTNacgtn]+$/) {
			print OUT $line."\n";
		} else {
			return "The file $file is not a FASTA-formatted file or contains a character other than A, C, G, T or N!\n";
		}
	}
	#################
	
	# close output file handler
	close OUT;
	# close input file handler
	close IN;
	
	return $output_dir."/".$basename.".tmp";
}
##############################

##############################
# split SAM-formatted files into mapped (flag 0 and 16) and unmapped (any other flag numbers) reads
# input: SAM-formatted files, directory for mapped reads, directory for unmapped reads, pattern for reference
# output: error message or 1
sub split_SAM {
	#################
	# save input parameters
	my @files = @{$_[0]};
	my $dir_mapped = $_[1];
	my $dir_unmapped = $_[2];
	my $pattern = $_[3];	
	#################
	
	# flag; generate list ids of all references
	my $first = 1;
	
	# open reference file handler
	open (REF, ">", $dir_mapped.$pattern."_list.txt") or return "Cannot create file $dir_mapped${pattern}_list.txt!\n";
	
	#################
	# for each SAM-file, save information about the alignment
	foreach my $file (@files) {
		# open input file handler
		open (IN, "<", $file) or return "Cannot open file $file!\n";
		
		my $basename = fileparse($file, @suffix);	
		# path to mapped file, extension '.map'
		my $map_file = $dir_mapped.$basename.".map";
		# path to unmapped file, extension '.unmap'
		my $unmap_file = $dir_unmapped.$basename.".rc";
		
		# open mapped file handler
		open (MAP, ">", $map_file) or return "Cannot create file $map_file!\n";
		# open unmapped file handler
		open (UNMAP, ">", $unmap_file) or return "Cannot create file $unmap_file!\n";

		while (<IN>) {
			my $line = $_;
			if ($line !~ /^\@/) {
				my @splitted_line = split(/\t/, $line);
				# usually a SAM fomat line starting with '@' has at least 10 columns
				if (scalar(@splitted_line) > 10) {
					# id of query
					my $read = $splitted_line[0];
					# mapping flag
					my $flag = $splitted_line[1];
					# reference id
					my $ref = $splitted_line[2];
					# start position of the query in the reference sequence
					my $start = $splitted_line[3];
					# CIGAR string
					my $cigar = $splitted_line[5];
					# query sequence
					my $seq = $splitted_line[9];
					
					# assumption: if the read id contains 'count', the file was in rc-format
					if ($read =~ /count/) {
						$read =~ s/\_/ /g;
					}
					
					#################
					# consider only those alignments with flag 0 or 16
					# flag = 0: alignment on forward strand
					# flag = 16: alignment on reverse strand
					if ($flag eq 0 || $flag eq 16 || $flag eq 256 || $flag eq 272) {
						my $mismatches = 0;
						my $strand = "+";
						
						#################
						# handle read mapped to reverse strand 
						if ($flag eq 16  || $flag eq 272) {
							$strand = "-";
							$seq = reverse($seq);
							$seq =~ tr/ACGTNacgtn/TGCANtgcan/;
						}
						#################
						
						#################
						# get mismatches
						my $mis_pos_str = "";
						if ($line =~ /MD:Z:(.*)?\t/) {
							my $md = $1;
							$mismatches = ($md =~ tr/[a-z|A-Z]//);
							my $tmp_pos = 0;
							if ($mismatches > 0 ) {
								while ($md =~ /((\d+)\^*(\D))/) {
									my $pattern = $1;
									my $pos = $2;
									my $nt = $3;
									$tmp_pos += $pos + 1;
									$mis_pos_str .= "," if (length($mis_pos_str) > 0);
									$mis_pos_str .= $pos.",";
									if ($pattern =~ /\^/) { # Bowtie 2 runs with the option: no gaps within first 32nts
										$mis_pos_str .= $tmp_pos."^".$nt;
										$md =~ s/$pos\^$nt//;
									} else {
										$mis_pos_str .= $tmp_pos."->".$nt;
										$md =~ s/$pattern//;
									}
								}
							}
							$mis_pos_str .= "," if (length($mis_pos_str) > 0);
							$mis_pos_str .= $md;
						}
						#################
						
						# write mapped reads
						print MAP ">".$read." ref=".$ref." strand=".$strand." start=".$start." mismatches=".$mismatches." cigar=".$cigar." align=".$mis_pos_str."\n";
						print MAP $seq."\n";
					} else {
						# write unmapped reads
						print UNMAP ">".$read."\n";
						print UNMAP $seq."\n";
					}
					#################
				}
			} elsif ($first && $line =~ /^\@SQ\sSN:(.*)\tLN:(\d+)/) {
				# save list of reference ids and length of the sequences
				print REF $1."\tlength=".$2."\n";
			}
		}
		# close mapped file handler
		close MAP;
		# close unmapped file handler
		close UNMAP;
		# close input file handler
		close IN;
		
		# generate list of references once
		$first = 0 if ($first);
	}
	#################
	
	# close reference file handler
	close REF;
	
	return 1;
}
##############################

#loaded ok				
1;