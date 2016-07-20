#!/usr/bin/perl
#03/29/2012
#Author: Petra Stepanowsky

package DiffExpr;

use warnings;
use strict;
use File::Basename;
use Statistics::R;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = qw(diffexpression);
@EXPORT_OK   = qw(diffexpression);
%EXPORT_TAGS = (DEFAULT  => [qw(&diffexpression)]);

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
# differential expression for 2 groups with at least 2 samples per group using DESeq
# input: table formatted file, path of output directory, hash with all group names and number of samples in group,
# groups used for differential expression, number of samples per group for differential expression
# output: error message or 1
sub diffexpression {
	#################
	# save input parameters
	my $file = $_[0];
	my $output_dir = $_[1];
	my %group_names = %{$_[2]};
	my @groups = @{$_[3]};
	my @num_samples = @{$_[4]};
	#################
	
	my @all_groups = keys %group_names;
	
	#################
	# there must be two groups with at least two samples per group
	if (scalar(@groups) == 2 && scalar(@num_samples) == 2 && $num_samples[0] >= 2 && $num_samples[1] >= 2) {
		# create output directory
		my $res = create_directory($output_dir); return $res if ($res ne 1);
	
		# create a communication bridge with R
		my $R = Statistics::R->new();
		
		# load DESeq library
		$R->run(q'library(DESeq)');

		
		#################
		# generate vector with group labels according to number of samples
		$R->send(qq'group1 <- rep("$groups[0]", $num_samples[0])');
		$R->send(qq'group2 <- rep("$groups[1]", $num_samples[1])');
		# if there are three groups, extract name of third group
		if (scalar(@all_groups) == 3) {
			my %used_groups = map{$_ => 1} @groups;
			my @left_group = grep(!defined $used_groups{$_}, @all_groups);
			if (scalar(@left_group) > 0) {
				$R->send(qq'group3 <- rep("$left_group[0]", $group_names{$left_group[0]})');
			}
		}
		#################
		
		$R->send(qq'output_dir <- "$output_dir"');
		$R->send(qq'aggregate_file <- "$file"');
		
		# initializes file names
		$R->run(q'unnorm_count_file <- paste(output_dir, "unnorm_count.txt", sep="")');
		$R->run(q'norm_count_file <- paste(output_dir, "norm_count.txt", sep="")');
		$R->run(q'result_file_g2_over_g1 <- paste(output_dir, "Group_", group2[1], "_over_", "Group_", group1[1], "_DESeq.txt", sep="")');
		$R->run(q'result_file_g1_over_g2 <- paste(output_dir, "Group_", group1[1], "_over_", "Group_", group2[1], "_DESeq.txt", sep="")');
		
		# read aggregate file
		$R->run(q'aggregate_count <- read.table(aggregate_file, header=TRUE)');
		$R->run(q'rownames(aggregate_count) <- aggregate_count[,1]');
		
		# sum rows
		$R->run(q'row_sum <- rowSums(aggregate_count[,2:ncol(aggregate_count)])');
		$R->run(q'col_sum <- colSums(aggregate_count[,2:ncol(aggregate_count)])');
		$R->run(q'tmp <- cbind(aggregate_count, row_sum)');
		$R->run(q'aggregate_count <- tmp[tmp[,ncol(tmp)]>0,1:ncol(tmp)-1]');
		$R->run(q'row_sum <- rowSums(aggregate_count[,2:ncol(aggregate_count)])');
		
		# unnormalized count file
		$R->run(q'unnorm_count<-cbind(aggregate_count, row_sum)');
		# sort by row sum
		$R->run(q'unnorm_count<-unnorm_count[order(unnorm_count[,dim(unnorm_count)[2]], decreasing=TRUE),]');
		# write to file
		$R->run(q'write.table(unnorm_count, unnorm_count_file, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)');
		
		# normalized count file
		$R->run(q'size_factors <-  col_sum/10000');
		$R->run(q'num_mirna <- dim(aggregate_count)[1]');
		$R->run(q'norm_count <- cbind(aggregate_count[1], round(t(t(aggregate_count[-1]) / size_factors), 7))');
		$R->run(q'row_sum <- rowSums(norm_count[-1])');
		$R->run(q'norm_count <- cbind(norm_count, row_sum)');
		$R->run(q'norm_count <- norm_count[order(norm_count[,ncol(norm_count)], decreasing=TRUE),]');
		$R->run(q'rownames(norm_count) <- aggregate_count[,1]');
		$R->run(q'colnames(norm_count) <- colnames(unnorm_count)');
		$R->run(q'write.table(norm_count, norm_count_file, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)');
		
		
		$R->run(q'conds <- c(group1, group2)');
		
		if (scalar(@all_groups) == 3) {
			
			# select columns of aggregate_count
			$R->run(q'index_group1 <- grep(paste(group1[1], "_", sep=""), colnames(aggregate_count))');
			$R->run(q'index_group2 <- grep(paste(group2[1], "_", sep=""), colnames(aggregate_count))');
			$R->run(q'index_group3 <- grep(paste(group3[1], "_", sep=""), colnames(aggregate_count))');
		} else {
			
			# select columns of aggregate_count
			$R->run(q'index_group1 <- grep(paste(group1[1], "_", sep=""), colnames(aggregate_count))');
			$R->run(q'index_group2 <- grep(paste(group2[1], "_", sep=""), colnames(aggregate_count))');
		}
		
		# indizes of group columns
		$R->run(q'index <- c(index_group1, index_group2)');

		$R->run(q'tmp <- aggregate_count[,index]');
		$R->run(q'row_sum <- rowSums(tmp[,2:ncol(tmp)])');
		$R->run(q'tmp <- cbind(tmp, row_sum)');
		$R->run(q'aggregate_count <- tmp[tmp[,ncol(tmp)]>0,2:ncol(tmp)-1]');
		
		$R->run(q'conds <- sort(conds)');
		$R->run(q'cds <- newCountDataSet(aggregate_count, conds)');
		$R->run(q'sizeFactors(cds) <- col_sum[index-1]/10000');
		
		# differential expression
		#$R->run(q'cds <- estimateVarianceFunctions(cds)');
                # the function estimateVarianceFunctions had been removed in DESeq 1.8.3 compiled for R 2.15.1
                #    use  estimateDispersions function instead as below 
                # $R->run(q'cds <- estimateDispersions(cds)');
	        ### With error catching
                $R->run(q'cds_out = try(estimateDispersions(cds),silent=TRUE)');
                $R->run(q'if (class(cds_out) == "try-error"){cds_out <- estimateDispersions(cds,fitType="local")}');
                $R->run(q'cds <- cds_out');		


		# g2 over g1
		$R->run(q'res <- nbinomTest(cds, group1[1], group2[1])[,1:8]');
		$R->run(q'res[,2:8] <- round(res[,2:8], 7)');
		$R->run(q'colnames(res)[3:5] <- c(paste("group_", group1[1], "_mean", sep=""), paste("group_", group2[1], "_mean", sep=""), paste("FC_", group2[1], "_over_", group1[1], sep=""))');
		$R->run(q'res[do.call(cbind,lapply(res,is.nan))] <- 0');
		$R->run(q'write.table(res, result_file_g2_over_g1, row.names=F, col.names=T, sep="\t", quote=F)');
		
		# g1 over g2
		$R->run(q'res <- nbinomTest(cds, group2[1], group1[1])[,1:8]');
		$R->run(q'res[,2:8] <- round(res[,2:8], 7)');
		$R->run(q'colnames(res)[3:5] <- c(paste("group_", group2[1], "_mean", sep=""), paste("group_", group1[1], "_mean", sep=""), paste("FC_", group1[1], "_over_", group2[1], sep=""))');
		$R->run(q'res[do.call(cbind,lapply(res,is.nan))] <- 0');
		$R->run(q'write.table(res, result_file_g1_over_g2, row.names=F, col.names=T, sep="\t", quote=F)');
				
		$R->stop();

		return 1;
	} else {
		return "Differential expression with less than 2 samples per group is not possible!\n";
	}
	#################
}
##############################

#loaded ok				
1;
