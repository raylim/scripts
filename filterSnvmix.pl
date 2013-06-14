#!/usr/bin/perl -w

use strict;

my $thresh = $ARGV[0];
my $minreads = $ARGV[1];

#chr10:100003467 C A C:1,A:1,0.7263945623,0.2604054213,0.0132000164,1
#chr10:171178    C       A       C:6,A:1,0.9912986338,0.0087013108,0.0000000553,1

while (<STDIN>) {
    chomp;
    my ($chrpos, $ref, $nonref, $notes) = split(/\s+/,$_);
    my ($refNum, $nonrefNum, $paa, $pab,$pbb,$call) = split (/,/,$notes);
    my ($R, $numR) = split(/:/, $refNum);
    my ($N, $numN) = split(/:/, $nonrefNum);
#    my ($numR,$R,$numN,$N,$paa,$pab,$pbb,$call) = split(/,/,$notes);
    if ((($numR+$numN)>$minreads) && (($pab+$pbb)>$thresh)) {
	print $_ ."\n";
    }
}

