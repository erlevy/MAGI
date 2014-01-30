#!/usr/bin/perl
#created 03/29/2012, last updated 08/10/2012
#Author: Petra Stepanowsky

package TargetPred;

use warnings;
use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = qw(targetprediction);
@EXPORT_OK   = qw(targetprediction);
%EXPORT_TAGS = (DEFAULT  => [qw(&targetprediction)]);

##############################
# predicts target sites for given microRNA reads using miRanda or CUDA-miRanda
# input: FASTA-formatted input file, FASTA-formatted reference file, score threshold, directory of executables (called CUDAMiranda or miranda)
# path and name of the output file, GPU enabled flag, reference batch size, number of available GPUs
# output: error code or 1
sub targetprediction {
	my $input_file = $_[0];
	my $reference = $_[1];
	my $score = $_[2];
	my $exe_dir = $_[3];
	my $output_file = $_[4];
	my $gpu = $_[5];
	my $rbs = $_[6];
	my $number_gpus = $_[7];
	my $max4gpu = 900;
	
	my $res = 0;
	if (-s $input_file) {
		if ($gpu eq "yes") {
			# CUDA-miRanda
			#$res = system("$exe_dir/CUDAMiranda $input_file $reference -sc $score -quiet -out $output_file -keyval -nGPUs $number_gpus -rbs $rbs -trim 30000");
			my $total = get_total_lines($input_file); 
			my $number_splitfiles = ($total - ($total % $max4gpu) ) / $max4gpu + 1;

			for(my $i=0; $i<$number_splitfiles; $i++) {
				$res = system("$exe_dir/CUDAMiranda $input_file.split.$i $reference -sc $score -quiet -out $output_file.split.$i -keyval -nGPUs $number_gpus -rbs $rbs -trim 30000");
			}

			$res = system("cat $output_file.split.* >> $output_file");

		} else {
			# miRanda
			$res = system("$exe_dir/miranda $input_file $reference -sc $score -quiet -out $output_file -keyval");
		}
	}
	
	$res == 0 ? return 1 : return $res;
}

# get the total line number of a file
sub get_total_lines{
	my $filename = $_[0];
	my $cnt = 0;

	open(FH, "<", $filename) or die "Cannot open a file $filename. $!";
	$cnt++ while <FH>;
	close FH;

	return $cnt; 
}

#loaded ok				
1;
