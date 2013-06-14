#!/usr/bin/env perl
# usage: perl getReference.pl {samtools} {bam-file}

use strict;

my $samtools = $ARGV[0];

my $bam = $ARGV[1];

my @header = `$samtools view -H $bam`;

my @SQ = grep /^\@SQ/, @header;

my $fasta_path = "/projects/pubseq/" . ($SQ[0] =~ m:downloads/(genomes[^\s]+):)[0];
print $fasta_path . "\n";

my @fasta_files = glob("$fasta_path/*.fasta $fasta_path/*.fa");

print $bam . " : " . $fasta_files[0] . "\n";
