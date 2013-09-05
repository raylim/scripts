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

    my $upstream = ($F{upstream_gene} eq $F{gene_name1})? 1 : 2;
    my $downstream = ($F{downstream_gene} eq $F{gene_name1})? 1 : 2;

    my $upstreamChr = $F{"gene_chromosome" . $upstream };
    my $downstreamChr = $F{"gene_chromosome" . $upstream };

    my $upstreamPosn = $F{"genomic_break_pos" . $upstream };
    my $downstreamPosn = $F{"genomic_break_pos" . $upstream };

    print join("\t", ($upstreamChr, $upstreamPosn, $downstreamChr, $downstreamPosn)) . "\n";
}
