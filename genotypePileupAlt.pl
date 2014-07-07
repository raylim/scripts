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
-o: output prefix for depth, alt depth, and af
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

my @names;
for my $bamFile (@bamFiles) {
    my $n = $bamFile;
    $n =~ s/.*\///;
    $n =~ s/\..*//;
    push @names, $n;
}

open DP, ">$opt{o}.dp.txt";
open AF, ">$opt{o}.af.txt";
open ADP, ">$opt{o}.altdp.txt";

print DP "CHROM\tPOS\tALT\t" . join("\t", @names) . "\n";
print ADP "CHROM\tPOS\tALT\t" . join("\t", @names) . "\n";
print AF "CHROM\tPOS\tALT\t" . join("\t", @names) . "\n";
for my $chrom (keys %{$chromPosAlts}) {
    for my $pos (keys %{$chromPosAlts->{$chrom}}) {
        my $region = "$chrom:$pos-$pos";
        for my $alt (@{$chromPosAlts->{$chrom}->{$pos}}) {
            print ADP "$chrom\t$pos\t$alt";
            print DP "$chrom\t$pos\t$alt";
            print AF "$chrom\t$pos\t$alt";
            for my $bamFile (@bamFiles) {
                my $nAlt = 0;
                my $dp = 0;
                my $af = 0;
                my @lines = `$samtools mpileup -r $region -f $opt{f} $bamFile 2> /dev/null`;
                if (@lines != 0) {
                    my $samOut = $lines[0];
                    chomp $samOut;
                    my @F = split /\t/, $samOut;
                    if (@F > 5) {
                        $dp = $F[3];
                        my $readBases = uc $F[4];
                        my $nAlt = $readBases =~ tr/$alt//;
                        $af = $nAlt / $dp if ($dp > 0);
                    }
                }
                print AF "\t$af";
                print DP "\t$dp";
                print ADP "\t$nAlt";
            }
            print AF "\n";
            print DP "\n";
            print ADP "\n";
        }
    }
}
close DP;
close ADP;
close AF;
