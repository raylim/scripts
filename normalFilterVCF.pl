#!/usr/bin/env perl
# filter tumor based on normal
# usage: normalFilterVCF.pl [tumor.vcf] [normal.vcf]

use strict;

if (@ARGV != 2) {
    print "Usage: normalFilterVCF.pl [tumor.vcf] [normal.vcf]\n" and exit(1);
}

my $tumorVCF = $ARGV[0];
my $normalVCF = $ARGV[1];

my $varPosn = {};
open IN, $normalVCF or die("Unable to open " . $normalVCF . "\n");
while (<IN>) {
	next if /^#/;
	my @F = split /\t/;
	my $chr = $F[0];
	my $posn = $F[1];
	$varPosn->{$chr} = {} unless exists $varPosn->{$chr};
	$varPosn->{$chr}{$posn} = 1;
}
close IN;

open IN, $tumorVCF or die("Unable to open " . $tumorVCF . "\n");
while (<IN>) {
	print if /^#/;
	my @F = split /\t/;
	my $chr = $F[0];
	my $posn = $F[1];
	print unless exists $varPosn->{$chr}{$posn};
}
close IN;
