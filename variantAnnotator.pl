#!/usr/bin/env perl # process snpEff output to pick the highest impact variant 
use strict;

my $chrPosn = {};
while (<>) {
    print $_ and next if /^#/;
    my @F = split /\t/;
    my $chr = $F[0];
    my $posn = $F[1];
    $chrPosn->{$chr} = {} unless exists $chrPosn->{$chr};
    $chrPosn->{$chr}{$posn} = [] unless exists $chrPosn->{$chr}{$posn};
    push @{$chrPosn->{$chr}{$posn}}, \@F;
}

my @highImpactEffects = qw/SPLICE_SITE_ACCEPTOR  SPLICE_SITE_DONOR START_LOST EXON_DELETED FRAME_SHIFT STOP_GAINED STOP_LOST/;
my @modImpactEffects = qw/NON_SYNONYMOUS_CODING CODON_CHANGE CODON_INSERTION CODON_CHANGE_PLUS_CODON_INSERTION CODON_DELETION CODON_CHANGE_PLUS_CODON_DELETION UTR_5_DELETED UTR_3_DELETED/;
my @lowImpactEffects = qw/SYNONYMOUS_START NON_SYNONYMOUS_START START_GAINED SYNONYMOUS_CODING SYNONYMOUS_STOP NON_SYNONYMOUS_STOP/;
my @modifiers = qw/NONE CHROMOSOME CUSTOM CDS GENE TRANSCRIPT EXON INTRON_CONSERVED UTR_5_PRIME UTR_3_PRIME DOWNSTREAM INTRAGENIC INTERGENIC INTERGENIC_CONSERVED UPSTREAM REGULATION INTRON/;

my %impactMap;
for my $effect (@highImpactEffects) {
    $impactMap{$effect} = 4;
}
for my $effect (@modImpactEffects) {
    $impactMap{$effect} = 3;
}
for my $effect (@lowImpactEffects) {
    $impactMap{$effect} = 2;
}
for my $effect (@modifiers) {
    $impactMap{$effect} = 1;
}

sub parseEffect {
    my $s = $_[0];
    $s =~ s/:.*//;
    return $s;
}

foreach my $chr (keys %{$chrPosn}) {
    foreach my $posn (keys %{$chrPosn->{$chr}}) {
        for my $F (@{$chrPosn->{$chr}{$posn}}) {
            print STDERR $F->[15] . " does not map" unless exists $impactMap{parseEffect($F->[15])};
        }
        my @sortedF = sort { $impactMap{parseEffect($a->[15])} <=> $impactMap{parseEffect($b->[15])} } @{$chrPosn->{$chr}{$posn}};
        #print $impactMap{parseEffect($_->[15])} . ":\t" . join("\t", @$_) for (@sortedF);
        my $topF = pop @sortedF;
        print join "\t", @$topF;
    }
}
