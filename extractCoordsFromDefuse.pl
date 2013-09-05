#!/usr/bin/env perl

use strict;
use warnings;

my @header;
while (my $line = <>) {
    chomp $line;

    if ($line =~ /^cluster_id/) {
        @header = split /\t/, $line;
        next;
    }

    my @arr = split /\t/, $line;
    my %F = map { $_ => shift @arr } @header;

    print "gene1: " . $F{gene_name1} . "\n";
    print "gene2: " . $F{gene_name2} . "\n";
    print "upstream gene: " . $F{upstream_gene} . "\n";
    print "downstream gene: " . $F{downstream_gene} . "\n";
    my $upstream = ($F{upstream_gene} eq $F{gene_name1})? 1 : 2;
    my $downstream = ($F{downstream_gene} eq $F{gene_name1})? 1 : 2;

    my $upstreamChr = $F{"gene_chromosome" . $upstream };
    my $downstreamChr = $F{"gene_chromosome" . $upstream };

    my $upstreamPosn = $F{"genomic_break_pos" . $upstream };
    my $downstreamPosn = $F{"genomic_break_pos" . $upstream };

    print join("\t", ($upstreamChr, $upstreamPosn, $downstreamChr, $downstreamPosn)) . "\n";
}
