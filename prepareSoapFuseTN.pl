#!/usr/bin/env perl
# prepare soapFuse file structure using samples.txt file
# print out soapfuse samples file

use strict;
use warnings;

use File::Path qw/make_path/;

use Getopt::Std;
my %opt;
getopts('h', \%opt);

my $usage = <<ENDL;
Usage: ./prepareSoapFuse.pl [samples_sets.txt]
ENDL

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}

HELP_MESSAGE if $opt{h};

sub getReadLength {
    my $fqFile = $_[0];
    open IN, "zcat $fqFile | head |" or die "Unable to open $fqFile\n";
    <IN>;
    return length(<IN>);
}


while (my $line = <>) {
    chomp $line;
    my @sampleSet = split / /, $line;
    # last sample is normal
    my $normalSample = $sampleSet[$#sampleSet];
    for my $tumorSample (@sampleSet[0..($#sampleSet-1)]) {
        my $tFq1 = "fastq/$tumorSample.1.fastq.gz";
        my $tFq2 = "fastq/$tumorSample.2.fastq.gz";
        my $nFq1 = "fastq/$normalSample.1.fastq.gz";
        my $nFq2 = "fastq/$normalSample.2.fastq.gz";
        die "Cannot find fastq files ($tFq1 and $tFq2)" unless (-e $tFq1 && -e $tFq1);
        die "Cannot find fastq files ($nFq1 and $nFq2)" unless (-e $nFq1 && -e $nFq1);
        my $tSampleDir = "soapfuse/$sample-T/$sample";
        my $nSampleDir = "soapfuse/$sample-N/$sample";
        make_path($tSampleDir);
        make_path($nSampleDir);
        system "ln -f $tFq1 $tSampleDir/${sample}_1.fastq.gz";
        system "ln -f $tFq2 $tSampleDir/${sample}_2.fastq.gz";
        system "ln -f $nFq1 $nSampleDir/${sample}_1.fastq.gz";
        system "ln -f $nFq2 $nSampleDir/${sample}_2.fastq.gz";
        my $tReadLength = &getReadLength($tFq1);
        my $nReadLength = &getReadLength($nFq1);
        print "$sample-T\t$sample\t$sample\t$tReadLength\n";
        print "$sample-N\t$sample\t$sample\t$nReadLength\n";
    }
}



