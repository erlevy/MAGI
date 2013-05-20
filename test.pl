#!/usr/bin/perl

my $input_file = "/storage/kawasakidisease/magirun/gpuyes/combined/novel/novel.fa";
my $total = get_total_lines($input_file);
print "Total number of lines is ".$total."\n";

my $number_splitfiles = ($total - ($total % 900) ) / 900 + 1;
print "The number of split files is ".$number_splitfiles."\n";




# get the total line number of a file
sub get_total_lines{
        my $filename = $_[0];
        my $cnt = 0;
        
        open(fh, $filename) or die "File open error. $!";
        $cnt++ while <fh>;
        close fh;

        return $cnt;
}

