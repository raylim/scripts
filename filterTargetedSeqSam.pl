#!/usr/bin/env perl
# filter G->T mutation artefacts from name sorted sam file

use strict;
use warnings;

use Getopt::Std;
my %opt;
getopts('h', \%opt);

my $usage = <<ENDL;
perl filterTargetedSeqSam.pl < [namesorted sam]
ENDL

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}

HELP_MESSAGE if $opt{h};

# print ref -> alt mismatches
sub parseMD {
    my ($seq, $md) = @_;
    my $origSeq = $seq;
    my $origMD = $md;
    #print "md: $md\n";
    my @refs;
    my @alts;
    do {
        $md =~ s/^([0-9]+)//;
        my $nmatch = $1;
        #print "nmatch: $nmatch\n";
        #kprint "seq: $seq\n";
        #if ($nmatch > length($seq)) {
       #    print "ALERT\n";
    #    print $origSeq . "\n";
#        print $origMD . "\n";
#        }
        $seq = substr $seq, $nmatch;
        if ($md =~ s/^([ATCGN])//) {
            my $ref = $1;
            my $alt = substr $seq, 0, 1;
            $seq = substr $seq, 1;
            #print "$ref $alt\n";
            push @refs, $ref;
            push @alts, $alt;
        }
        $md =~ s/^\^[NATCG]+//;
    } while ($md ne "");
    (\@refs, \@alts);
}

sub countGT {
    my ($refs, $alts) = @_;
    my $count = 0;
    for my $i (0..$#$refs) {
        if ($refs->[$i] eq "G" && $alts->[$i] eq "T") {
            $count++;
        }
    }
    $count;
}


my @gTab;
while (defined(my $p1 = <>) and defined(my $p2 = <>)) {
    chomp $p1;
    chomp $p2;
    my @F1 = split /\t/, $p1;
    my @F2 = split /\t/, $p2;
    my @optF1 = @F1[11..$#F1];
    my @optF2 = @F2[11..$#F2];
    my @mdF1s = grep /^MD:Z:/, @optF1;
    my @mdF2s = grep /^MD:Z:/, @optF2;
    next if (@mdF1s < 1 || @mdF2s < 1);
    my $mdF1 = $mdF1s[0];
    my $mdF2 = $mdF2s[0];
    $mdF1 =~ s/^MD:Z://; 
    $mdF2 =~ s/^MD:Z://; 
    my ($refs1, $alts1) = parseMD($F1[9], $mdF1);
    my ($refs2, $alts2) = parseMD($F2[9], $mdF2);
    my $count1 = countGT($refs1, $alts1);
    my $count2 = countGT($refs2, $alts2);
    #print $count1 + $count2 . "\t$count1\t$count2\n";
    $gTab[$count1 + $count2]++;
}
my $total;
for my $i (0..$#gTab) {
    next unless defined $gTab[$i];
    $total += $gTab[$i];
}
for my $i (0..$#gTab) {
    next unless defined $gTab[$i];
    printf "%i\t%i\t%.2f\n", $i, $gTab[$i], $gTab[$i] / $total;
}
print "Total\t$total\n";
