#!/usr/bin/env perl
# filter PE fastq files based on number of Ns

use strict;
use warnings;

use Bio::SeqIO;
use Getopt::Std;

my %opt;
getopts('p:1:2:', \%opt);

my $usage = <<ENDL;
Usage: filterFastq.pl -p [max percent fraction of read] -1 output.1.fastq[.gz] -2 output.2.fastq[.gz] input.1.fastq[.gz] input.2.fastq[.gz]
-p [1-100]: maximium percentage fraction of the read that can be 'N'
-1 [file]: output file 1
-2 [file]: output file 2
ENDL

sub HELP_MESSAGE {
    print STDERR $usage;
    exit(1);
}

if (@ARGV != 2) {
    print STDERR "Need two fastq files\n" and HELP_MESSAGE();
}
unless ($opt{1} && $opt{2}) {
    print STDERR "Need two output files\n" and HELP_MESSAGE();
}
unless ($opt{p}) {
    print STDERR "Need minimum fracton of read\n" and HELP_MESSAGE();
}




my $fq1 = $ARGV[0];
$fq1 = "gunzip -c $fq1 |" if ($fq1 =~ /\.gz$/);
my $fq2 = $ARGV[1];
$fq2 = "gunzip -c $fq2 |" if ($fq2 =~ /\.gz$/);

my $outfq1 = ($opt{1} =~ /\.gz$/)? "| gzip -c > $opt{1} " : "> $opt{1}";
my $outfq2 = ($opt{2} =~ /\.gz$/)? "| gzip -c > $opt{2} " : "> $opt{2}";

my $seqin1 = Bio::SeqIO->new(
    -file => $fq1,
    -format => 'fastq');
my $seqin2 = Bio::SeqIO->new(
    -file => $fq2,
    -format => 'fastq');

my $seqout1 = Bio::SeqIO->new(
    -file => $outfq1,
    -format => 'fastq');
my $seqout2 = Bio::SeqIO->new(
    -file => $outfq2,
    -format => 'fastq');

my $i = 0;
my $j = 0;
while ((my $seqobj1 = $seqin1->next_seq) && (my $seqobj2 = $seqin2->next_seq)) {
    my $threshold = length($seqobj1->seq) * ($opt{p} / 100);
    my $n1 = ($seqobj1->seq =~ tr/N//);
    my $n2 = ($seqobj2->seq =~ tr/N//);

    unless ($n1 > $threshold || $n2 > $threshold) {
        $seqout1->write_seq($seqobj1);
        $seqout2->write_seq($seqobj2);
        $j++;
    }
    $i++;
    printf STDERR ("Processed %i sequences: wrote %i (%i filtered)\n", $i, $j, $i-$j) if $i % 100000 == 0;
}

printf STDERR "Finished: Processed %i sequences; wrote %i; filtered %i\n", $i, $j, $i - $j;

