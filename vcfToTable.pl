#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Std;
my %opt;
getopts('hg:', \%opt);

#my $usage = <<ENDL;

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}

HELP_MESSAGE if $opt{h};

my $cmd = $extractFields . 

