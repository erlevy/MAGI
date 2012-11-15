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
	
	my $res = 0;
	if (-s $input_file) {
		if ($gpu eq "yes") {
			# CUDA-miRanda
			$res = system("$exe_dir/CUDAMiranda $input_file $reference -sc $score -quiet -out $output_file -keyval -nGPUs $number_gpus -rbs $rbs -trim 30000");
		} else {
			# miRanda
			$res = system("$exe_dir/miranda $input_file $reference -sc $score -quiet -out $output_file -keyval");
		}
	}
	
	$res == 0 ? return 1 : return $res;
}

#loaded ok				
1;