# 5/2/2014
# Author: Eric Levy
# Test command:
#perl -MExportResults -e'ExportResults::exportresults("/storage/magi/acuvir/");'

#!/usr/bin/perl
package ExportResults;

use warnings;
use strict;
use File::Basename;
use File::Copy::Recursive qw(dircopy);
use File::Copy;
use Cwd;

sub exportresults {
	#################
	# save input parameters
	# an existing output directory
	my $dir = $_[0];

	#################
	# specify output directories
	my $dir_rc = $dir."rc/";
	my $dir_visualization = $dir."visualization/";
	my $dir_visualization_downsized = $dir_visualization."preprocessing/";
	my $dir_visualization_mapping = $dir_visualization."alignment/";
	my $dir_visualization_known = $dir_visualization."known/";
	my $dir_visualization_novel = $dir_visualization."novel/";
	my $dir_visualization_target = $dir_visualization."target/";

	my $dir_output = $dir."output/";
	my $dir_output_rc = $dir_output."rc/";
	my $dir_output_downsized = $dir_output."preprocessing/";
	my $dir_output_mapping = $dir_output."alignment/";
	my $dir_output_known = $dir_output."known/";
	my $dir_output_novel = $dir_output."novel/";
	my $dir_output_target = $dir_output."target/";
	#################
	create_directory($dir_output);
	
	if (-d $dir_rc) {
		dircopy($dir_rc, $dir_output_rc) or return "Cannot copy quality files!\n";
	}

	if (-d $dir_visualization_downsized) {
		dircopy($dir_visualization_downsized, $dir_output_downsized) or return "Cannot copy quality files!\n";
	}

	if (-d $dir_visualization_mapping) {
		dircopy($dir_visualization_mapping, $dir_output_mapping) or return "Cannot copy quality files!\n";
	}

	if (-d $dir_visualization_known) {
		dircopy($dir_visualization_known, $dir_output_known) or return "Cannot copy quality files!\n";
	}
	if (-d $dir_visualization_novel) {
		dircopy($dir_visualization_novel, $dir_output_novel) or return "Cannot copy quality files!\n";
	}
	if (-d $dir_visualization_target) {
		create_directory($dir_output_target);
		copy($dir_visualization_target."novel_vs_genes.txt", $dir_output_target."novel_vs_genes.txt") or return "Cannot copy novel.fa!\n";
	}
	if (-d $dir."combined/") {
		create_directory($dir_output."combined/");
		dircopy($dir."combined/alignment/", $dir_output."combined/alignment/") or return "Cannot copy quality files!\n";
		dircopy($dir."combined/mapped/", $dir_output."combined/mapped/") or return "Cannot copy quality files!\n";
		if (-e $dir."combined/combined_aggregate.bound"){	
			copy($dir."combined/combined_aggregate.bound", $dir_output."combined/combined_aggregate.bound") or return "Cannot copy novel.fa!\n";
		}	
		if (-d $dir."combined/novel/") {
			create_directory($dir_output."combined/novel/");
			copy($dir."combined/novel/novel.fa", $dir_output."combined/novel/novel.fa") or return "Cannot copy novel.fa!\n";
			copy($dir."combined/novel/precursors.str", $dir_output."combined/novel/precursors.str") or return "Cannot copy novel.fa!\n";
			copy($dir."combined/novel/probs_precursors.txt", $dir_output."combined/novel/probs_precursors.txt") or return "Cannot copy novel.fa!\n";
		}
	}
	my $old_directory = cwd;
	chdir("$dir") or die "can't change directory";
	my $res = system("zip -r output.zip"." output/");
	chdir($old_directory);

	return 1;
}

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

1;
