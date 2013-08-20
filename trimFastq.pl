#!/usr/bin/env perl
# trim fastq file

use strict;
use warnings;

use Getopt::Std;

my %opt;
getopts('hl:', \%opt);

my $usage = <<ENDL;
Usage: filterFastq.pl -l [read length] 
-h: this help message
-l [integer]: length
ENDL

sub HELP_MESSAGE {
    print STDERR $usage;
    exit(1);
}

print "Missing read-trim length\n" and HELP_MESSAGE() unless $opt{l};

my $i = 0;
while (<STDIN>) {
    chomp;
    if ($i % 2 == 0) {
        print;
    } else {
        print substr($_, 0, $opt{l});
    }
    print "\n";
    $i++;
}
