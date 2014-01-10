#!/usr/bin/env perl
# zero base quality scores for G->T and C->A mutation artefacts from sorted sam file

use strict;
use warnings;

use Getopt::Std;
my %opt;
getopts('h', \%opt);

my $usage = <<ENDL;
perl fixTargetedSeqSam.pl < [namesorted sam]
ENDL

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}

HELP_MESSAGE if $opt{h};

my %gtFilter = (
    "xT" => 1,
    "xC" => 1);

my %caFilter = (
    "xT" => 1,
    "xC" => 1);

#my %gtFilter = (
#"TxG" => 1,
#"TxA" => 1,
#"GxA" => 1,
#"CxG" => 1,
#"CxA" => 1,
#"GxA" => 1,
#"AxA" => 1
#);

#my %gtFilter = (
#"TxGT" => 1,
#"TxAT" => 1,
#"TxAT" => 1,
#"TxA" => 1,
#"CxG" => 1,
#"CxA" => 1,
#"AxA" => 1,
#"AxG" => 1
#);

#my %caFilter = (
#"TxT" => 1,
#"TxG" => 1,
#"TxC" => 1,
#"TxA" => 1,
#"CxG" => 1,
#"CxA" => 1
#);



# extract ref -> alt mismatches and contexts using MD 
sub parseMD {
    my ($seq, $md) = @_;
    my $origSeq = $seq;
    #my $origMD = $md;
    #print "md: $md\n";
    my @refs;
    my @alts;
    my @contexts;
    my @posns;
    my $i = 0;
    do {
        $md =~ s/^([0-9]+)//;
        my $nmatch = $1;
        $i += $nmatch;
        $seq = substr $seq, $nmatch;
        if ($md =~ s/^([ATCGN])//) {
            my $ref = $1;
            my $alt = substr $seq, 0, 1;
            $i++;
            push @posns, $i;
            # pad sequence
            my $ss = " " . $origSeq . " " ;
            substr($ss, $i, 1) = "x";
            my $context = substr $ss, $i, 2;
            $seq = substr $seq, 1;
            #print "$ref $alt\n";
            push @refs, $ref;
            push @alts, $alt;
            push @contexts, $context;
        }
        $md =~ s/^\^[NATCG]+//;
    } while ($md ne "");
    (\@refs, \@alts, \@posns, \@contexts);
}

while (defined(my $line = <>)) {
    print $line and next if $line =~ /^\@/;
    chomp $line;
    my @F = split /\t/, $line;
    my $seq = $F[9];
    my $bq = $F[10];

    my @optF = @F[11..$#F];
    my @mdFs = grep /^MD:Z:/, @optF;
    next if (@mdFs < 1);
    my $mdF = $mdFs[0];
    $mdF =~ s/^MD:Z://; 
    my ($refs, $alts, $posns, $contexts) = parseMD($seq, $mdF);
    for my $i (0..$#$refs) {
        if ($refs->[$i] eq "G" && $alts->[$i] eq "T" && exists $gtFilter{$contexts->[$i]}) {
            #print substr($seq, $posns->[$i] - 1, 1) . "\n";
            substr($bq, $posns->[$i] - 1, 1) = "!";
            #    print $contexts->[$i] . "\n";
        }
        if ($refs->[$i] eq "C" && $alts->[$i] eq "A" && exists $caFilter{$contexts->[$i]}) {
            substr($bq, $posns->[$i] - 1, 1) = "!";
        }
    }
    $F[10] = $bq;

    print join("\t", @F) . "\n";
}
