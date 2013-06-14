#!/usr/bin/perl
# Convert mutect tables to vcf files

use strict;

# print header
my $now = localtime;
my $header = <<END;
##fileformat=VCFv4.1
##fileDate=$now
##FILTER=<ID=covered,Description="site powered to detect a mutation (80% power for a 0.3 allelic fraction mutation)">
##INFO=<ID=POW,Number=A,Type=Float,Description="tumor_power * normal_power">
##INFO=<ID=TPOW,Number=A,Type=Float,Description="given the tumor sequencing depth, what is the power to detect a mutation at 0.3 allelic fraction">
##INFO=<ID=NPOW,Number=A,Type=Float,Description="given the normal sequencing depth, what power did we have to detect (and reject) this as a germline variant">
##INFO=<ID=DP,Number=1,Type=Integer,Description="total tumor and normal read depth which come from paired reads">
##INFO=<ID=IMP,Number=1,Type=Integer,Description="number of reads which have abnormal pairing (orientation and distance)">
##INFO=<ID=Q0,Number=1,Type=Integer,Description="total number of mapping quality zero reads in the tumor and normal at this locus">
##INFO=<ID=LOGLIK,Number=A,Type=Float,Description="Log of (likelihood tumor event is real / likelihood event is sequencing error )">
##INFO=<ID=AF,Number=A,Type=Float,Description="allelic fraction of this candidated based on read counts">
##INFO=<ID=CF,Number=A,Type=Float,Description="estimate of contamination fraction used (supplied or defaulted)">
##INFO=<ID=CONTAMINANTLOGLIK,Number=A,Type=Float,Description="log likelihood of ( event is contamination / event is sequencing error )">
##FORMAT=<ID=REFDP,Number=1,Type=Integer,Description="count of reference alleles">
##FORMAT=<ID=ALTDP,Number=1,Type=Integer,Description="count of alternate alleles">
##FORMAT=<ID=REFSUM,Number=1,Type=Integer,Description="sum of quality scores of reference alleles">
##FORMAT=<ID=ALTSUM,Number=1,Type=Integer,Description="sum of quality scores of alternate alleles">
##FORMAT=<ID=INS,Number=1,Type=Integer,Description="number of insertion events at this locus">
##FORMAT=<ID=DEL,Number=1,Type=Integer,Description="number of deletion events at this locus">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=NLOGLIK,Number=A,Type=Float,Description="log likelihood of ( normal being reference / normal being altered )">
END

print $header;


my @header = split /\t/,<>;

while (<>) {



my %info;
my %tumorFormat;
my %normalFormat;

