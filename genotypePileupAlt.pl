#!/usr/bin/env perl
# run samtools mpileup on positions to determine allelic frequency and depth
# input format: chrom, pos, allele

use strict;
use warnings;

use FIle::Temp qw/ tempfile tempdir /;

use Getopt::Std;

my %opt;
getopts('hs:f:o:', \%opt);

my $usage = <<ENDL;
Usage: genotypePileupAlt.pl -o [output prefix] -s [samtools] -f [ref fasta file] [chrom-pos-alt file] [bam files]
-h: this help message
-s: samtools binary to use
-f: reference file
-o: output prefix for depth and af
ENDL

print "Missing ref fasta\n" and HELP_MESSAGE() unless $opt{f};
print "Missing output prefix\n" and HELP_MESSAGE() unless $opt{o};

my $samtools = 'samtools';
$samtools = $opt{s} if $opt{s};

my $cpaFile = shift @ARGV;
my @bamFiles = @ARGV;

my $tmpBed = File::Temp->new(SUFFIX => '.bed');
open IN, $cpaFile;
my @alts;
my @chromPosAlts;
my $i = 0;
while (<IN>) {
    chomp;
    my @F = split /\t/;
    print $tmpBed $F[0] . "\t" . $F[1] - 1;
    $chromPosAlts[$i] = join "\t", ($F[0], $F[1], $F[2]);
    $alts[$i] = $F[2];
    $i++;
}

my %af;
my %dp;
my @names;
for my $bamFile (@bamFiles) {
    my $tmp = File::Temp->new();
    system("$samtools mpileup -f $opt{f} $bamFile > $tmp->filename");
    my $n = $bamFile;
    $n =~ s/.*\///;
    $n =~ s/\..*//;
    push @names, $n;

    my $i = 0;
    while (<$tmp>) {
        chomp;
        my @F = split /\t/;
        my $cov = $F[3];
        my $readBases = uc $F[4];
        my $alt = uc $alts[$i++];
        my $nAlt = $readBases =~ tr/$alt//;
        my $dp{$n}[$i] = $cov;
        my $af{$n}[$i] = $nAlt / $cov;
        $i++;
    }
}

open DP, ">$opt{o}.dp.txt";
open AF, ">$opt{o}.af.txt";

print "CHROM\tPOS\tALT\t" . join("\t", @names) . "\n";
for my $n (@names) {
    for my $i (0 .. $#alts) {
        print DP $chromPosAlt[$i] . "\t" . $dp{$n}[$i] . "\n";
        print AF $chromPosAlt[$i] . "\t" . $af{$n}[$i] . "\n";
    }
}

close DP;
close AF;

