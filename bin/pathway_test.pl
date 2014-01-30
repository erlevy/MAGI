#!/usr/bin/perl

use strict;
use warnings;
use PathwayEnrichment;

my $output_dir = "/storage/magi/6a51a039e5363dba7f4802fc88678f0a947389ce/";
my $path_data = "/home/triton/pipelineTyler/data/pathway/";
my $path_bin = "/home/triton/pipelineTyler/bin/";

##### perform pathway enrichment #####
	print "..... perform pathway enrichment: ";
	PathwayEnrichment::pathwayenrichment($output_dir, $path_data, $path_bin);
	print "done\n";
#####
