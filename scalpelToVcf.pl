#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Std;
my %opt;
getopts('hf:n:t:', \%opt);

my $usage = <<ENDL;
perl scalpelToVcf.pl -f [ref fasta] -t [tumor sample name] -n [normal sample name]
ENDL

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}


HELP_MESSAGE if $opt{h};

my $now = localtime;
my $normal = $opt{n};
my $tumor = $opt{t};

open REF, $opt{f} . ".fai" or die "Unable to open reference index file\n";

sub getRefSeq {
    my ($chrom, $pos) = @_;
    open(SAMTOOLS, "samtools faidx $opt{f} $chrom:$pos-$pos |");
    <SAMTOOLS>;
    my $seq = <SAMTOOLS>;
    chomp $seq;
    close(SAMTOOLS);
    $seq;
}

my @contigs = <REF>;

my $vcfHeader = <<ENDL;
##fileformat=VCFv4.1
##fileDate=$now
##INFO=<ID=INDELTYPE,Number=1,Type=String,Description="insertion or deletion">
##INFO=<ID=SIZE,Number=1,Type=Integer,Description="size of variant">
##INFO=<ID=AVGKCOV,Number=1,Type=Float,Description="Average k-mer coverage">
##INFO=<ID=MINKCOV,Number=1,Type=Integer,Description="Minimum k-mer coverage">
##INFO=<ID=CHI2SCORE,Number=1,Type=Float,Description="Chi-2 score">
##INFO=<ID=COVRATIO,Number=1,Type=Float,Description="Coverage ratio">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=1,Type=String,Description="K-mer allelic depth">
##FORMAT=<ID=AF,Number=1,Type=String,Description="K-mer variant allelic frequency">
ENDL
##PEDIGREE=<Derived=$tumor,Original=$normal>
for my $contig (@contigs) {
    my @F = split /\t/, $contig;
    $vcfHeader .= "##contig=<ID=$F[0],length=$F[1]>\n";
}
$vcfHeader .= "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$tumor\t$normal\n";
print $vcfHeader;

my $format = "GT:AD:AF";

my @header;
while (<>) {
    chomp;
    if (/^#chr/) {
        $_ =~ s/^#//;
        @header = split/\t/;
        next;
    }
    next if /^#/;
    print "$_\n";
    my @F = split /\t/;
    my %F = map { $_ => shift @F } @header;
    my $chrom = $F{chr};
    my $pos = $F{start};
    my $ref = $F{ref};
    my $alt = $F{obs};
    if ($ref eq "-" || $alt eq "-") {
        $pos--;
        my $seq = getRefSeq($chrom, $pos);
        $alt = $seq . $alt;
        $ref = $seq . $ref;
        $alt =~ s/-//;
        $ref =~ s/-//;
    }
    my ($normalR, $tumorR, $normalA, $tumorA) = $F{bestState} =~ /(\d)(\d)\/(\d)(\d)/;
    my $normalGT = "0" x $normalR . "1" x $normalA;
    $normalGT = join "/", (split '', $normalGT);
    my $tumorGT = "0" x $tumorR . "1" x $tumorA;
    $tumorGT = join "/", (split '', $tumorGT);

    my ($normalRcov, $tumorRcov, $normalAcov, $tumorAcov) = $F{covState} =~ /^(\d+) (\d+)\/(\d+) (\d+)/;
    my $normalAD = "$normalRcov,$normalAcov";
    my $tumorAD = "$tumorRcov,$tumorAcov";
    my $normalAF = $normalAcov / ($normalRcov + $normalAcov);
    my $tumorAF = $tumorAcov / ($tumorRcov + $tumorAcov);

    my $tumorFormat = sprintf "%s:%s:%.2f", $tumorGT, $tumorAD, $tumorAF;
    my $normalFormat = sprintf "%s:%s:%.2f", $normalGT, $normalAD, $normalAF;

    my $info = "SIZE=$F{size};INDELTYPE=$F{type};AVGKCOV=$F{avgKcov};MINKCOV=$F{minKcov};CHI2SCORE=$F{chi2score};COVRATIO=$F{covRatio}";
    print "$chrom\t$pos\t.\t$ref\t$alt\t.\tPASS\t$info\t$format\t$tumorFormat\t$normalFormat\n";
}
