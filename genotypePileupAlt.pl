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

my $tmpBed = File::Temp->new();
open IN, $cpaFile;
my $chromPosAlts = {};
while (<IN>) {
    chomp;
    my ($chrom, $pos, $alt) = split "\t";
    $chromPosAlts->{$chrom} = {} unless exists $chromPosAlts->{$chrom};
    $chromPosAlts->{$chrom}->{$pos} = [] unless exists $chromPosAlts->{$chrom}->{$pos};
    push @{$chromPosAlts->{$chrom}->{$pos}}, uc($alt);
}
close IN;

my $adp = {};
my $dp = {};
my @names;
for my $bamFile (@bamFiles) {
    my $n = $bamFile;
    $n =~ s/.*\///;
    $n =~ s/\..*//;
    push @names, $n;
    $dp->{$n} = {};
    $adp->{$n} = {};

    print "$samtools view -L $cpaFile -f $opt{f} -b $bamFile | $samtools mpileup -d 20000 -BQ0 -l $cpaFile -f $opt{f} - 2> /dev/null\n";
    my @lines = `$samtools view -L $cpaFile -f $opt{f} -b $bamFile | $samtools mpileup -d 20000 -BQ0 -l $cpaFile -f $opt{f} - 2> /dev/null`;
    for my $line (@lines) {
        chomp $line;
        my @F = split /\t/, $line;
        next if (scalar(@F) != 6);
        my $chrom = $F[0];
        my $pos = $F[1];
        for my $alt (@{$chromPosAlts->{$chrom}->{$pos}}) {
            my $cov = length($F[4]);
            my $readBases = uc $F[4];
            my $nAlt = $readBases =~ tr/$alt//;
            $dp->{$n}->{$chrom}->{$pos}->{$alt} = $cov;
            $adp->{$n}->{$chrom}->{$pos}->{$alt} = $nAlt;
        }
    }
}

open DP, ">$opt{o}.dp.txt";
open ADP, ">$opt{o}.adp.txt";
open AF, ">$opt{o}.af.txt";

print DP "CHROM\tPOS\tALT\t" . join("\t", @names) . "\n";
print ADP "CHROM\tPOS\tALT\t" . join("\t", @names) . "\n";
print AF "CHROM\tPOS\tALT\t" . join("\t", @names) . "\n";
for my $chrom (keys %{$chromPosAlts}) {
    for my $pos (keys %{$chromPosAlts->{$chrom}}) {
        for my $alt (@{$chromPosAlts->{$chrom}->{$pos}}) {
            print DP "$chrom\t$pos\t$alt";
            print ADP "$chrom\t$pos\t$alt";
            print AF "$chrom\t$pos\t$alt";
            for my $n (@names) {
                my $cov = 0;
                my $nAlt = 0;
                $cov = $dp->{$n}->{$chrom}->{$pos}->{$alt} if (exists $dp->{$n}->{$chrom}->{$pos}->{$alt});
                $nAlt = $adp->{$n}->{$chrom}->{$pos}->{$alt} if (exists $adp->{$n}->{$chrom}->{$pos}->{$alt});
                print DP "\t$cov";
                print ADP "\t$nAlt";
                print AF "\t" . (($cov > 0)? $nAlt / $cov : 0);
            }
            print DP "\n";
            print ADP "\n";
            print AF "\n";
        }
    }
}
close DP;
close ADP;
close AF;
