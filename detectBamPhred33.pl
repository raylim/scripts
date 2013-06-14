#!/usr/bin/perl
# usage : samtools view [bam file] | ./detectBamPhred33.pl 
# returns: exit code 0 if bam is phred 33 score format

use strict;

while (my $line = <>) {
    my @F = split /\t/, $line;
    if (@F[10] =~ /[!"#$%&'()*+,\-.\/0123456789:]/) {
        # phred 33
        print "Phred 33 detected\n";
        exit(0);
    } elsif (@F[10] =~ /[JKLMNOPQRSTUVWXYZ[\\\]^_`abcdefghijklmnopqrstuvwxyz{|}~]/) {
        # phred 64
        print "Phred 64 detected\n";
        exit(1);
    }
}

exit(1);
