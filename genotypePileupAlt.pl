#!/usr/bin/env perl
# run samtools mpileup on positions to determine allelic frequency and depth
# input format: chrom, pos, allele

use strict;
use warnings;

use File::Temp qw/ tempfile tempdir /;

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
my %chromPosAlts;
my @cpaKeys;
while (<IN>) {
    chomp;
    my @F = split "\t";
    print $tmpBed $F[0] . "\t" . ($F[1] - 1) . "\n";
    my $cp = join ":", ($F[0], $F[1]);
    push @{$chromPosAlts{$cp}}, uc($F[2]);
    push @cpaKeys, join(":", @F);
}
close IN;

my %af;
my %dp;
my @names;
for my $bamFile (@bamFiles) {
    my $tmp = File::Temp->new();
    my $sys = "$samtools mpileup -l $tmpBed -f $opt{f} $bamFile > $tmp";
    print $sys . "\n";
    system("$sys") == 0 or die "samtools failed : $?"; 
    my $n = $bamFile;
    $n =~ s/.*\///;
    $n =~ s/\..*//;
    push @names, $n;

    $dp{$n} = {};
    $af{$n} = {};
    my $i = 0;
    while (<$tmp>) {
        print;
        chomp;
        my @F = split /\t/;
        my $cov = $F[3];
        my $readBases = uc $F[4];
        for my $alt (@{$chromPosAlts{$cp}}) {
            my $nAlt = $readBases =~ tr/$alt//;
            $dp{$n}->{$cpa}, $cov;
            $af{$n}->{$cpa}, $nAlt / $cov;
        }
    }
}

open DP, ">$opt{o}.dp.txt";
open AF, ">$opt{o}.af.txt";

print "chromPosAlt\t" . join("\t", @names) . "\n";
for my $cpa (@cpaKeys) {
    print DP "$cpa";
    print AF "$cpa\t";
    for my $n (@names) {
        my $d = (exists $dp{$n}->{$cpa})? $dp{$n}->${cpa} : 0;
        my $a = (exists $af{$n}->{$cpa})? $af{$n}->${cpa} : 0;
        print DP "\t$d";
        print AF "\t$a";
    }
    print DP "\n";
    print AF "\n";
}

close DP;
close AF;

