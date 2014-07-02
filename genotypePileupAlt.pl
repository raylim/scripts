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
my $chromPosAlts = {};
while (<IN>) {
    chomp;
    my ($chrom, $pos, $alt) = split "\t";
    $chromPosAlts->{$chrom} = {} unless exists $chromPosAlts->{$chrom};
    $chromPosAlts->{$chrom}->{$pos} = [] unless exists $chromPosAlts->{$alt};
    push @{$chromPosAlts->{$chrom}->{$pos}}, uc($alt);
}
close IN;

my $af = {};
my $dp = {};
my @names;
for my $bamFile (@bamFiles) {
    my $n = $bamFile;
    $n =~ s/.*\///;
    $n =~ s/\..*//;
    push @names, $n;
    $dp->{$n} = {};
    $af->{$n} = {};
    for my $chrom (keys %{$chromPosAlts}) {
        $dp->{$n}->{$chrom} = {};
        $af->{$n}->{$chrom} = {};
        for my $pos (keys %{$chromPosAlts->{$chrom}}) {
            $dp->{$n}->{$chrom}->{$pos} = {};
            $af->{$n}->{$chrom}->{$pos} = {};
            my $region = "$chrom:$pos-$pos";
            my @lines = `$samtools mpileup -r $region -f $opt{f} $bamFile 2> /dev/null`;
            if (@lines == 0) {
                # no sam output
                for my $alt (@{$chromPosAlts->{$chrom}->{$pos}}) {
                    $dp->{$n}->{$chrom}->{$pos}->{$alt} = 0;
                    $af->{$n}->{$chrom}->{$pos}->{$alt} = 0;
                }
            } else {
                my $samOut = $lines[0];
                chomp $samOut;
                my @F = split /\t/, $samOut;
                my $cov = $F[3];
                my $readBases = uc $F[4];
                for my $alt (@{$chromPosAlts->{$chrom}->{$pos}}) {
                    my $nAlt = $readBases =~ tr/$alt//;
                    $dp->{$n}->{$chrom}->{$pos}->{$alt} = $cov;
                    $af->{$n}->{$chrom}->{$pos}->{$alt} = $nAlt / $cov;
                }
            }
        }
    }
}

open DP, ">$opt{o}.dp.txt";
open AF, ">$opt{o}.af.txt";

print DP "CHROM\tPOS\tALT\t" . join("\t", @names) . "\n";
print AF "CHROM\tPOS\tALT\t" . join("\t", @names) . "\n";
for my $chrom (keys %{$chromPosAlts}) {
    for my $pos (keys %{$chromPosAlts->{$chrom}}) {
        for my $alt (@{$chromPosAlts->{$chrom}->{$pos}}) {
            print DP "$chrom\t$pos\t$alt";
            print AF "$chrom\t$pos\t$alt";
            for my $n (@names) {
                print DP "\t" . $dp->{$n}->{$chrom}->{$pos}->{$alt};
                print AF "\t" . $af->{$n}->{$chrom}->{$pos}->{$alt};
            }
            print DP "\n";
            print AF "\n";
        }
    }
}
close DP;
close AF;
