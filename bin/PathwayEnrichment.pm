
# Test command
#perl -MPathwayEnrichment -e'PathwayEnrichment::pathwayenrichment("/storage/magi/6a51a039e5363dba7f4802fc88678f0a947389ce/", "../data/pathway/", "./");'

#!/usr/bin/perl
package PathwayEnrichment;

use warnings;
#use strict;
use Statistics::R;


sub pathwayenrichment {

	print "Beginning Pathway Enrichment\n";
	
	my $output_dir = $_[0];
	my $pathway_dir = $_[0]."pathway/";
	my $pathway_genes_dir = $pathway_dir."genes/";
	my $data_dir = $_[1];
	my $bin_dir = $_[2];

	print $output_dir."\n";
	print $data_dir."\n";
	print $bin_dir."\n";

	my $novel_target_input = $output_dir."visualization/target/novel_vs_genes.txt";
	open my $input_mir, '<', $novel_target_input;

	my $precursor_input = $output_dir."combined/novel/probs_precursors.txt";
	open my $input_prec, '<', $precursor_input;

	my $genes_ref_input = $data_dir."microtcds_targets_intersect_unique.txt";
	open my $genes_ref, '<', $genes_ref_input;

	unless(-e $pathway_dir or mkdir $pathway_dir) 
	{
			die "Unable to create $pathway_dir\n";
	}

	unless(-e $pathway_genes_dir or mkdir $pathway_genes_dir) 
	{
			die "Unable to create $pathway_genes_dir\n";
	}

	my $novel_target_output = $pathway_dir."novel_vs_genes_filtered.txt";
	open my $target_output, '>', $novel_target_output;

	my $precursor_output = $pathway_dir."probs_precursors_filtered.txt";
	open my $prec_output, '>', $precursor_output;

	# OUTLINE
	# input: probs_precursors.txt
	# actions: sort, take top (20?)
	# output: list of top 20 as "novel#"

	# input: novel_vs_genes.txt
	# actions: take just entries for top 20 novel, sort, take top (10?) 
	# output: list of top 10 targets for each novel

	my @genes = <$genes_ref>;
	s/\s+\z// for @genes;

	my $probs_threshold = 0.95;
	my $target_threshold = 500;

	my $prec_num = 0;
	my $target_num = 0;
	my @target_mirs = ();
	#my @targets = ();
	while ( <$input_prec> )
	{
		chomp;
		my (@columns) = split /\s+/, $_;
		if ($columns[1] > $probs_threshold) 
		{
			my $out_line = $columns[0] . "	" . $columns[1] . "	\n";
			print $prec_output $out_line;
			$target_num += 1;
			push(@target_mirs, substr $columns[0], 1);
		}
		$prec_num += 1;
	}
	print "Filtered from $prec_num to $target_num novel target_mirs.\n";

	my $n1 = 0;
	my $n2 = 0;
	#my @targets = ();
	while ( <$input_mir> )
	{
		chomp;
		my (@columns) = split /\s+/, $_;
		next if ($columns[0] eq "miRNA");
	#	print "$columns[2]\n";
	#	print "@mir\n";
		if (substr($columns[0], -1, 1) eq '*')
		{
			$mir = $columns[0];
			chop $mir;
			$mir = substr $mir, 5;
		}
		else
		{
			$mir = $columns[0];
			$mir = substr $mir, 5;
		}
		if (($mir ~~ @target_mirs) && $columns[2] > $target_threshold && ($columns[1] ~~ @genes))
		{
	#		print "@columns\n";
	#		push(@targets,\@columns);
	#		unless ($columns[1] ~~ @targets)
	#		{
	#			push(@targets, $columns[1]);
	#		}
			my $out_line = $columns[0] . "	" . $columns[1] . "	" . $columns[2] . "	\n";
			print $target_output $out_line;
			$n2 += 1;
		}
		$n1 += 1;
	}
	print "Filtered from $n1 to $n2 miR-target pairs.\n";

	my $known_diff_input = $output_dir."known/diff_expr/Group_g1_over_Group_g2_DESeq.txt";
	open my $known_diff, '<', $known_diff_input;

	my $novel_diff_input = $output_dir."novel/diff_expr/Group_g1_over_Group_g2_DESeq.txt";
	open my $novel_diff, '<', $novel_diff_input; 

	my $known_top_output = $pathway_dir."known_top.txt";
	open my $known_top_genes, '>', $known_top_output;

	my $novel_top_output = $pathway_dir."novel_top.txt";
	open my $novel_top_genes, '>', $novel_top_output;

	# OUTLINE
	# input: diffexp for known and novel
	# actions: take all below p-value threshold (0.05?)
	# output: list of all known and novel miR with diffexp below threshold

	my $p_threshold = 0.05;
	my $line = "\_";

	my $mir_input_num = 0;
	my $n4 = 0;
	my @input_mirs = ();
	#my @targets = ();
	my @novel = ();
	while ( <$novel_diff> )
	{
		chomp;
		my (@columns) = split /\s+/, $_;
		next if ($columns[0] eq "id");
		if ($columns[7] < $p_threshold)
		{
			my $out = $columns[0] . "\n";
			#print $novel_top_genes $out_line;
			#$n4 += 1;
			#push(@input_mirs, substr $columns[0], 1);
			print $novel_top_genes $out;
			push(@novel, $columns[0]);
		}
		$mir_input_num += 1;
	}

	my @known = ();
	while ( <$known_diff> )
	{
		chomp;
		my (@columns) = split /\s+/, $_;
		next if ($columns[0] eq "id");
		if ($columns[7] < $p_threshold)
		{
			my $cut = length($columns[0]) - index($columns[0],$line);
			my $out = substr($columns[0], 0, -$cut)  . "\n";
			#print $novel_top_genes $out_line;
			#$n4 += 1;
			#push(@input_mirs, substr $columns[0], 1);
			print $known_top_genes $out;
			push(@known, substr($columns[0], 0, -$cut));
		}
		$mir_input_num += 1;
	}

	### miR List Inputs
	### Creates @mir of all miR in input list
	my $input1 = $pathway_dir."known_top.txt";
	open my $input_known, '<', $input1;

	my $input2 = $pathway_dir."novel_top.txt";
	open my $input_novel, '<', $input2;

	# Process miR lists input into arrays of miR names
	#my @known = ();
	#while ( <$input_known>)
	#{
	#	chomp;
	#	push(@known, $_);
	#	print "$_\n";
	#}

	#my @novel = ();
	#while ( <$input_novel>)
	#{
	#	chomp;
	#	push(@novel, $_);
	#}

	#print "@mir\n";

	### Test output file
	### Currently writes KEGG pathway gene table results
	my $output = $pathway_dir."pathway_tables.txt";
	open my $out, '>', $output;

	### microT Genes Input
	### Creates @targets
	### All target genes found in all known miR from input list
	my $file2 = $data_dir."microtcds_targets_intersect.txt";
	open my $microT, '<', $file2;

	# Process microT input into array of miR-target pairs
	# Only if the known miR is in the input list
	my @targets = ();
	while ( <$microT> )
	{
		chomp;
		my (@columns) = split /\s+/, $_;
	#	print "$columns[0]\n";
	#	print "@mir\n";
		if ($columns[0] ~~ @known)
		{
			unless ($columns[1] ~~ @targets)
			{
				push(@targets, $columns[1]);
			}
		}
	}

	#print "@targets\n";
	print "Number of gene targets from known: " . scalar (@targets) . "\n";

	### Novel Genes Input
	### Adds to @targets
	### All target genes found in all novel miR from input list
	my $file3 = $pathway_dir."novel_vs_genes_filtered.txt";
	open my $novelT, '<', $file3;

	# Process microT input into array of miR-target pairs
	# Only if the novel miR is in the input list
	while ( <$novelT> )
	{
		chomp;
		my (@columns) = split /\s+/, $_;
	#	print "$columns[0]\n";
	#	print "@mir\n";
		if ($columns[0] ~~ @novel)
		{
			unless ($columns[1] ~~ @targets)
			{
				push(@targets, $columns[1]);
			}
		}
	}

	#print "@targets\n";
	print "Number of gene targets from known and novel: " .scalar (@targets) . "\n";

	### KEGG Genes Input
	### Creates @kegg and @kegg_genes
	### List of all keg pathways, and list of references to corresponding genes in pathway
	my $file = $data_dir."kegg_intersect_ids_microt.txt";
	open my $in, '<', $file;

	# Separate arrays for kegg pathway names and list of references to arrays of the genes in that pathway
	my @kegg = ();
	my @kegg_id = ();
	my @kegg_genes = ();
	while (<$in>)
	{
		chomp;
		my (@columns) = split /\s+/, $_;
#		print "@columns\n";
		push(@kegg,$columns[0]);
		shift(@columns);
		push(@kegg_id,$columns[0]);
		shift(@columns);
		push(@kegg_genes,\@columns);
	}

#	foreach (@kegg)
#	{
#		print "$_\n";
#	}

#	foreach (@kegg_id)
#	{
#		print "$_\n";
#	}
#
#	foreach (@kegg_genes)
#	{
#		print "@$_\n";
#	}


	### Create count tables for all the KEGG pathways
	my @kegg_tables = ();
	my @kegg_targets = ();

	#print "@{$kegg_genes[0]}\n";

	for ( $i = 0; $i < scalar (@kegg); $i++ )
	{
		my $kegg_targets_output = $pathway_genes_dir."$kegg[$i].txt";
		open my $targets_out, '>', $kegg_targets_output;
		my @current_table = (0,0,0,0);
		my @current_kegg_genes = @{$kegg_genes[$i]};
	#	print "@current_kegg_genes\n";
		foreach (@targets)
		{
	#		print "@{$kegg_genes[$i]}\n";
	#		print "@current_table\n";
			if ($_ ~~ @current_kegg_genes)
			{
				$current_table[0]++;
				push(@kegg_targets,$_);
				print $targets_out "$_\n";
			}
		}
		$current_table[2] = scalar (@current_kegg_genes) - $current_table[0];
		push(@kegg_tables, \@current_table);
	}

	my $target = 0;
	my $non_target = 0;

	foreach (@kegg_tables)
	{
		my $current_table = $_;
		#print "@$current_table[0]\n";
		$target += @$current_table[0];
		$non_target += @$current_table[2];
	}

	foreach (@kegg_tables)
	{
		my $current_table = $_;
		@$current_table[1] = scalar (@targets) - @$current_table[0];
		@$current_table[3] = 4448 - scalar (@targets) - @$current_table[2];
	}

	my $n5 = 0;
	foreach (@kegg_tables)
	{
		my $current_table = $_;
		if (@$current_table[0] != 0)
		{
	#		print "@$_\n";
			$n5++; 
		}
	}
	print "Number of non-zero pathways: $n5\n";


	#my @values = ();

	#foreach (@kegg_tables)
	#{
	#	my $current_table = $_;
	#	my @numbers = (@$current_table[0], @$current_table[1], @$current_table[2], @$current_table[3]);
	#	push(@values, @$current_table[0]/@$current_table[1]);
	#}

	#@values = sort(@values);

	#foreach (@values)
	#{
	#	print "$_\n";
	#}

	#%values = ();
	#for ( $i = 0; $i < scalar (@kegg); $i++ )
	#{
	#	my $curr_pathway = $kegg[$i];
	#	my $curr_value = $values[$i];
	#	$values{$curr_pathway} = $curr_value;
	#}

	#@keys = sort {$values{$b} <=> $values{$a}} keys %values;

	#@key_slice = @keys[0..9];
	#foreach $key (@key_slice)
	#{
	#	print "$key: $values{$key}\n";
	#}

	#print "$target\n";
	#print "$non_target\n";
	for ( $i = 0; $i < scalar (@kegg); $i++ )
	{
		my $current_table = $kegg_tables[$i];
		my $table_out = @$current_table[0] . "	" . @$current_table[1] . "	" . @$current_table[2] . "	" . @$current_table[3];
		print $out $kegg[$i] . "	";
		print $out $kegg_id[$i] . "	";	
		print $out $table_out;
		print $out "\n";
	}

	my $path_table_file = $pathway_dir."pathway_tables.txt";
	my $path_out_file = $pathway_dir."all_pathways_out.txt";
	my $fishers_function_path = $bin_dir."/fishers_function.R";
	my $R = Statistics::R->new();
	$R->run(qq(b <- read.table("$path_table_file",stringsAsFactors=F)));
	$R->run(qq(source("$fishers_function_path")));
	$R->run(qq(write.table(output, "$path_out_file", sep = "	", quote = F, row.names = F, col.names = F)));

}

1;
