#!/usr/bin/perl -w 

use strict;
use File::Basename;
use DBI;

if (@ARGV != 2) {
    print "\tUsage:\n";
    print "\t\t./generateArtifactFinder <snvmix_novel_codon_annot file> <full path to output file>\n";
    exit;
}

#CONSTANT
my $EXTENSION = ".txt";
my $JUNCTIONSPLICE = "/share/data/rgoya/scripts/getIntronZone-beta.pl";

#command line args
my $nonsynfile = $ARGV[0];
my $artifactfile = $ARGV[1];

#files generated
my $rootpath = dirname($nonsynfile);
my $basename = basename($nonsynfile, $EXTENSION);
my $genenamefile = $rootpath . "/" . $basename . "_gene.txt";
my $tmpcheckfile = $rootpath . "/" . $basename . "_tmpCheck.txt";
my $tmpresultfile = $rootpath . "/" . $basename . "_tmpResult.txt";

#commands to run
my $juncArtCheck = "$JUNCTIONSPLICE -l $tmpcheckfile > $tmpresultfile";

# open up database to get gene name from ensemble id
#open a database handle

my $dbh = DBI->connect('DBI:mysql:ovmut_51;host=shannon.cluster.bccrc.ca;port=3306', 'ovmut_rw', 'o9v8m7')
    or die "Couldn't connect to the database: " . DBI->errstr;

# reusable prepared statement handles                                                                                                                                                                                                
my $gene_handle = $dbh->prepare("SELECT gene_name FROM Gene WHERE ensg_id=?");

# create temperory file for junction splicing check                                                                                                                                                                                  
open (INFO, $nonsynfile);
open (FILE, ">$tmpcheckfile");
open (FILE2, ">$genenamefile");

print FILE "gene_name\tchr\tstart\n";
while (my $line = <INFO>) {
    chomp($line);
    my @content = split(/\s+/, $line);
    my ($chromosome, $position) = split(/:/, $content[0]);
    $chromosome =~ s/chr//;

    $gene_handle->execute($content[1]);
    my $gene_result = $gene_handle->fetchrow();

    print FILE join("\t", $gene_result, $chromosome, $position) . "\n";

    print FILE2 join(" ", $line, $gene_result) . "\n";
}


close(INFO);
close(FILE);
close(FILE2);

# get junction splicing result                                                                                                                                                                                                       
#print $juncArtCheck . "\n";
`$juncArtCheck`;

# join junction splicing result with annot file                                                                                                                                                                                      
open (RESULT, $tmpresultfile);

my %splice_result = ();

print "$tmpresultfile\n";

while (my $newline = <RESULT>) {
    if ($newline =~ m/\bConnecting\b|^\#\#SUMMARY$|^\s*$/) {
        next;
    }
    chomp ($newline);
    my ($gene, $chr, $pos, $result) = split(/\t/, $newline);
    my $key = join ("\t", $gene, $chr, $pos);
    $splice_result{$key} = $result;
}

close (RESULT);

#print "size of hash: " . keys (%splice_result) . "\n";

open (MERGE, $genenamefile);
open (FINAL, ">$artifactfile");

while (my $read_content = <MERGE>) {
    chomp($read_content);
    my @info = split(/\s+/, $read_content);
    my ($mut_chr, $mut_pos)  = split(/:/, $info[0]);
    $mut_chr =~ s/chr//;
    my $mut_gene = pop(@info);

    my $mut_key = join("\t", $mut_gene, $mut_chr, $mut_pos);

    print FINAL join(" ", @info, $splice_result{$mut_key})  . "\n";
}

close (MERGE);
close (FINAL);

`rm $genenamefile $tmpcheckfile $tmpresultfile`;
