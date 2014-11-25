#!/usr/bin/env perl
# fix the museq vcf so that it conforms to spec

use strict;

use POSIX;
use feature qw(switch);

use Getopt::Std;

my %opt;
getopts('hR:', \%opt);

my $help = <<ENDL;
Usage: fixMuseqVCF.pl -R [reference.fasta] [vcf]
-R: reference.fasta with .dict
-h: this help message
ENDL

sub HELP_MESSAGE {
    print STDERR $help;
    exit(1);
}

HELP_MESSAGE() if $opt{h};
print STDERR "Missing reference.fasta\n" and HELP_MESSAGE() unless $opt{R};
print STDERR "Need museq vcf output to fix\n" and HELP_MESSAGE() unless @ARGV == 1;

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
##INFO=<ID=AD,Number=A,Type=Integer,Description="Allelic depths for the ref and alt">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=PR,Number=1,Type=Float,Description="Somatic probability">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=TriNTContext,Number=1,Type=String,Description="Tri-nucleotide context">
ENDL
for my $contig (@contigs) {
    $header .= "##contig=<ID=$contig,length=$contigLength{$contig}";
    $header .= ",assembly=$assembly" if $assembly;
    $header .= ">\n";
}
$header .= "##reference=file://$refFasta\n";

my $e = 0.0000001;
my $formatString = "GT:DP:AD";


open IN, $ARGV[0] or die "Unable to open $ARGV[0]\n";
my @lines = <IN>;
my @headerLines = grep /^#/, @lines;
my $tumorSampleFile;
my $normalSampleFile;
for my $line (@headerLines) {
    if ($line =~ /^##tumour=(.*)$/) {
        $tumorSampleFile = $1;
    }
    if ($line =~ /^##normal=(.*)$/) {
        $normalSampleFile = $1;
    }
}

my $normalSample = $normalSampleFile;
$normalSample =~ s/.*\///;
$normalSample =~ s/\..*//;
my $tumorSample = $tumorSampleFile;
$tumorSample =~ s/.*\///;
$tumorSample =~ s/\..*//;

#$header .= "##PEDIGREE=<Derived=$tumorSample,Original=$normalSample>\n";
$header .= "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$tumorSample\t$normalSample\n";
print STDOUT $header;

my @headerLines = grep(/^[^#]/, @lines);
foreach (@headerLines) {
    chomp;
    my @F = split /\t/;
    my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) = @F;
    next unless $ref =~ m/[ATCG]/;
    next unless $alt =~ m/[ATCG]/;
    my %infoMap;
    for (split /;/, $info) {
        my ($key, $value) = split /=/;
        $infoMap{$key} = $value;
    }
    my $tumorNRef = $infoMap{"TR"};
    my $tumorNAlt = $infoMap{"TA"};
    my $normalNRef = $infoMap{"NR"};
    my $normalNAlt = $infoMap{"NA"};
    my $nref = $tumorNRef + $normalNRef;
    my $nalt = $tumorNAlt + $normalNAlt;

    my $pr = $infoMap{"PR"};

    my $dp = $nref + $nalt;
    my $normalDP = $normalNAlt + $normalNRef;
    my $tumorDP = $tumorNAlt + $tumorNRef;

    my $ad = "$nref,$nalt";
    my $normalAD = "$normalNRef,$normalNAlt";
    my $tumorAD = "$tumorNRef,$tumorNAlt";

    $qual = floor($qual);
    #$qual = floor(-10 * log($pr) / log(10));

    my $normalGT = "0/0";
    my $tumorGT = ($tumorNRef == 0)? "1/1" : "0/1";

    $info = "DP=$dp;AD=$ad;PR=$pr";

    # GT:DP:AD
    my $tumorFormat = "$tumorGT:$tumorDP:$tumorAD";
    my $normalFormat = "$normalGT:$normalDP:$normalAD";

    print "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$formatString\t$tumorFormat\t$normalFormat\n";
}
close IN;

exit(0);
