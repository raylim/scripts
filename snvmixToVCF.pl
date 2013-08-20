#!/usr/bin/env perl
# convert snvmix output to vcf

use strict;

use POSIX;
use feature qw(switch);

use Getopt::Std;

my %opt;
getopts('d:hR:', \%opt);

my $help = <<ENDL;
Usage: snvmixToVCF.pl -R [reference.fasta] -d [depth filter] [table]
-R: reference.fasta with .dict
-d: non-ref depth filter threshold
-h: this help message
ENDL

sub HELP_MESSAGE {
    print STDERR $help;
    exit(1);
}

HELP_MESSAGE() if $opt{h};
print STDERR "Missing reference.fasta\n" and HELP_MESSAGE() unless $opt{R};
print STDERR "Need snvmix output to convert\n" and HELP_MESSAGE() unless @ARGV == 1;

my $refFasta = $opt{R};
my $refDict = $refFasta;
$refDict =~ s/\.[^.]+$/.dict/;
print STDERR "Missing reference.dict\n" and HELP_MESSAGE() unless (-e $refDict);

# read dictionary file
my @contigs;
my %contigLength;
open IN, $refDict or die("Unable to open $refDict\n");
my @dict = <IN>;
@dict = grep /^\@SQ/, @dict;
for my $line (@dict) {
    $line =~ m/SN:(.+)\tLN:(.+)\tUR:/;
    push @contigs, $1;
    $contigLength{$1} = $2;
}

my $assembly;
$assembly = "b37" if ("MT" ~~ @contigs);
$assembly = "hg19" if ("chrM" ~~ @contigs);

my $header = "##fileformat=VCFv4.1\n";
$header .= "##FILTER=<ID=AD,Description=\"allelic depth < $opt{d}\">\n" if $opt{d}; 
$header .= <<ENDL;
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt">
##FORMAT=<ID=DP,Number=.,Type=Integer,Description="Total depth of ref and alt">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
ENDL
for my $contig (@contigs) {
    $header .= "##contig=<ID=$contig,length=$contigLength{$contig}";
    $header .= ",assembly=$assembly" if $assembly;
    $header .= ">\n";
}
$header .= "##reference=file://$refFasta\n";

my $e = 0.0000001;
my $formatString = "GT:DP:AD:PL";

my $sampleName = $ARGV[0];
$sampleName =~ s/.*\///;
$sampleName =~ s/\..*//;
$header .= "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sampleName\n";

print STDOUT $header;
open IN, $ARGV[0] or die "Unable to open $ARGV[0]\n";
while (<IN>) {
    chomp;
    my @F = split /\t/;
    my ($chrpos, $ref, $alt, $csf) = @F[0..3];
    next unless $ref =~ m/[ATCG]/;
    next unless $alt =~ m/[ATCG]/;
    my ($chr, $pos) = $chrpos =~ m/(.+):([^\t]+)/;
    my ($nref, $nalt, $pAA, $pAB, $pBB, $maxP) = split /,/, $csf;
    $pAA += $e;
    $pAB += $e;
    $pBB += $e;
    $nref =~ s/^.://;
    $nalt =~ s/^.://;
    my $filter = ".";
    if ($opt{d}) {
        $filter = ($nalt > $opt{d})? "PASS" : "AD";
    }
    my $dp = $nref + $nalt;
    my $ad = "$nref,$nalt";
    my $qual = floor(-10 * log($pAA) / log(10));
    my $gt;
    given($maxP) {
        when(0) { $gt = "0/0"; }
        when(1) { $gt = "0/1"; }
        default { $gt = "1/1"; }
    }
    my $pl = sprintf "%.0f,%.0f,%.0f", -10 * log($pAA) / log(10), -10 * log($pAB) / log(10), -10 * log($pBB) / log(10);
    my $info = "DP=$dp";
    my $format = "$gt:$dp:$ad:$pl";
    print "$chr\t$pos\t.\t$ref\t$alt\t$qual\t$filter\t$info\t$formatString\t$format\n";
}
close IN;

exit(0);
