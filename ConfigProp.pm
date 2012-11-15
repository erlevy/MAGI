#!/usr/bin/perl
# 12/6/2011
# Author: Petra Stepanowsky

package ConfigProp;

use warnings;
use strict;
use Cwd 'abs_path';
use File::Basename;

our $VERSION = 1.00;

my $path = dirname(abs_path($0));

##############################
# returns a property by key
# input: key
# output: value or empty string
sub get_property {
	my $self = shift;
	if (keys %{$self->{properties}} == 0) {
		$self->init;
	}
	my $key = shift;
	if (exists $self->{properties}->{$key}) {
		return $self->{properties}->{$key};
	} else {
		return "";
	}
}
##############################

##############################
# initializes properties
# input: object of class ConfigProp
sub init {	
	my $self = shift;
	open (CONFIG, "<", $self->{config}) or die "Cannot open file $self->{config}!\n";
	while (<CONFIG>) {
		my $line = $_;
		if ($line !~ /^#/) {
			$line =~ s/\s//g;
			my @split = split(/=/, $line);
			if ($#split == 1) {
				chomp $split[1];
				$self->{properties}->{$split[0]} = $split[1];
			}
		}
	}
	close CONFIG;
}
##############################

##############################
# constructor
# input: configuration file
# output: object of class ConfigProp
sub new {	
	my($class, $config) = @_;
	my $self = {config => $config};
	bless($self, $class);
	$self->init;
	return $self;
}
##############################

#loaded ok				
1;