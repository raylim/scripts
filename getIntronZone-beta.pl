#!/usr/bin/perl
#
# Rodrigo Goya, 2009
# 
# Get information on the localization of a specific base pair in regards
# to the exons in a gene.
#
# Gets the information either gene by gene through the -g flag
# Gets the information in a list via tab-separated files using one of two formats:
#	gene_name	chr	pos
#	chr	pos
# The corresponding header HAS to be present
#
# It will get the information from ensembl, change information below to connect
# a local copy of ensembl-core database
#
# For each position it will report different character codes, indicating the position
# of the base-pair in regards to the exons in the gene
#
#	*	chr:pos was found inside an exon far from the junction sites
#	<	chr:pos was found inside an exon close to the 5' junction site
#	>	chr:pos was found inside an exon close to the 3' junction site
#	5	chr:pos was found outside an exon upstream, close to the 5' exon start
#	3	chr:pos was found outside an exon downstreaml close to the 3' exon end
#
# The "closeness" to the junction site:
#	- on the exon side, determined by the -e flag, default 10 bp
# 	- on the intron side, determined by the -i flag, default 5 bp
#
# In other 



use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;

# Try to connect, get it over with if we can't
print "Connecting to Ensembl core...\n";
#my $dbCore = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
#			-user   => "ensembl51_r",
#			-dbname => "homo_sapiens_core_65_37",
#			-host   => "shannon",
#			-pass   => "e4n3s2",			     
#			-driver => 'mysql',
#			-port   => 3306
#                );
# Change data base connection information here
# Beware of the port, ensembldb.ensembl.org uses 5306
# but many local installations use 3306
my $dbCore = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
			-user   => "anonymous",
			-dbname => "homo_sapiens_core_65_37",
			-host   => "ensembldb.ensembl.org",
			-pass   => "",			     
			-driver => 'mysql',
			-port   => 5306
               );

use Getopt::Std;
my $opt_string = 'hg:l:p:o:e:i:v';
my %opt;
getopts( "$opt_string", \%opt ) or usage();
usage() if $opt{h};

my %targets;
if($opt{g} && !$opt{l}) {
	my $geneID = $opt{g} if $opt{g} or usage();
	$geneID =~ tr/[a-z]/[A-Z]/;
	my $pos = 0;
	$pos = $opt{p} if $opt{p};
	push @{$targets{$geneID}{'pos'}},$pos;
} elsif($opt{l} && !$opt{g}) {
	open(INPUT,"<$opt{l}") || die("ERROR: cannot open file $opt{l} for reading\n");
	$_ = <INPUT>;
	if(m/^gene_name\s+chr\s+start\W/) {
		my %tmpPos;
		while(<INPUT>) {
			my ($geneID, $chr, $pos, $str, @tmp) = split(/\s+/, $_);
			$tmpPos{$geneID}{$pos} = 1;
			$tmpPos{$geneID}{'chr'} = $chr;
		}
		foreach my $geneID (keys %tmpPos) {
			$targets{$geneID}{'chr'} = $tmpPos{$geneID}{'chr'};
			foreach my $pos (keys %{$tmpPos{$geneID}}) {
				push @{$targets{$geneID}{'pos'}}, $pos unless $pos eq 'chr';
			}
		}
	} elsif(m/^chr\s+start\W/ || m/^chr\s+pos\W/) {
		my %tmpPos;
		my $slice_adaptor = $dbCore->get_SliceAdaptor;
		while(<INPUT>) {
			my ($chr, $pos, @tmp) = split(/\s+/, $_);
		 	my $slice = $slice_adaptor->fetch_by_region( 'chromosome', "$chr", $pos, $pos);
			my @genes = @{ $slice->get_all_Genes() };
			#foreach my $gene (@genes) { print "DEBUG BIOTYPE ".$gene->external_name."\t".$gene->biotype."\n"; }
			if(scalar(@genes) > 1) {
				my @selectGeneID;
				# Found more than two genes, both have CDS at that region?
				# Testing, check exons... this is inefficient, I know
				foreach my $gene (@genes) {
					foreach my $transcript ( @{$gene->get_all_Transcripts()} ) {
                        			foreach my $exon (@{ $transcript->get_all_Exons }) {
                        			        my $eStart = $exon->start();
                        			        my $eEnd = $exon->end();
	#print "DEBUG: ".$gene->slice->seq_region_name." ".$exon->stable_id." ".$exon->start()." <=0 && ".$exon->end()." >= 0\n";
							if( $exon->start() <= 0 && $exon->end() >= 0 ) {
								if( !exists {map { $_ => 1 } @selectGeneID}->{$gene->external_name} ) {
									push @selectGeneID, $gene->external_name;
								}
							}
						}
					}
				}
				if( scalar(@selectGeneID) > 1) {
					my $dump;
					foreach my $gene (@genes) {
						$dump .= "".$gene->slice->seq_region_name."\t".$gene->stable_id."\t".$gene->external_name."\n";
					}
					print("DEBUG: more than one gene found, dump:\n$dump\t-> picking ".$selectGeneID[0]."\n");
					
				}
				$tmpPos{$selectGeneID[0]}{$pos} = 1;
				$tmpPos{$selectGeneID[0]}{'chr'} = $chr;
			} elsif(scalar(@genes) == 1) {
				$tmpPos{$genes[0]->external_name}{$pos} = 1;
				$tmpPos{$genes[0]->external_name}{'chr'} = $chr;
			} else {
				print("Could not find genes at $chr:$pos\n");
			}
		}
		foreach my $geneID (keys %tmpPos) {
			$targets{$geneID}{'chr'} = $tmpPos{$geneID}{'chr'};
			foreach my $pos (keys %{$tmpPos{$geneID}}) {
				push @{$targets{$geneID}{'pos'}}, $pos unless $pos eq 'chr';
			}
		}
	} else {
		close(INPUT);
		die("ERROR: file $opt{l} unrecognizable\n");
	}
	close(INPUT);
} elsif($opt{g} && $opt{l}) {
	die("ERROR: cannot simultaneously give a gene and a list\n");
} else {
	usage();
}
my $fileOut;
$fileOut = $opt{o} if $opt{o};
# exon threshold
my $eThr = 10;
$eThr = $opt{e} if $opt{e};
# intron treshold
my $iThr = 5;
$iThr = $opt{i} if $opt{i};

my $gene_adaptor = $dbCore->get_GeneAdaptor();
my %summary;
foreach my $geneID (keys %targets) {
	my @genes = @{$gene_adaptor->fetch_all_by_external_name($geneID)};

	if($opt{v}) {
		print "".scalar(@genes)." matches found\n";
	}
	foreach my $gene (@genes) {
		
		if($opt{v}) {
			print "$geneID"; 
			if($opt{p} || $opt{l}) { print " (pos ".join(",",@{$targets{$geneID}{'pos'}}).")"; }
			print " ".get_string($gene,"s")."\n";
		}
		foreach my $transcript ( @{$gene->get_all_Transcripts()} ) {
			if($opt{v}) {
				print "\tTRANSCRIPT ".$transcript->stable_id()."\n";
			}
			foreach my $exon (@{ $transcript->get_all_Exons }) {
				my $eStart = $exon->start();
				my $eEnd = $exon->end();
				my $posStr;
if($eStart > $eEnd) { 
	print "Start: $eStart > End: $eEnd !!!\n";
}
				if($opt{p} || $opt{l}) {
					foreach my $pos (@{$targets{$geneID}{'pos'}}) {
#print "DEBUG: ".$exon->stable_id." $pos >= ".$exon->start()." && $pos <= ".$exon->end()."\n";
				   		if( ($pos >= $eStart) && ($pos <= $eEnd) ) {
							if( ($pos > $eStart+$eThr) && ($pos < $eEnd-$eThr) ) {
								$posStr .= "*";
								$summary{$geneID}{$pos}{'*'}++;
							} else {
								if($pos >= $eEnd-$eThr) {
									$posStr .= ">";
									$summary{$geneID}{$pos}{'>'}++;
								}
								if($pos <= $eStart+$eThr) {
									$posStr .= "<";
									$summary{$geneID}{$pos}{'<'}++;
								}
							}
						} elsif( ($pos > $eStart-$iThr) && ($pos < $eStart)  ) {
							$posStr .= "5";
							$summary{$geneID}{$pos}{'5'}++;
						} elsif( ($pos < $eEnd+$iThr) && ($pos > $eEnd)  ) {
							$posStr .= "3";
							$summary{$geneID}{$pos}{'3'}++;
						}
					}
				}
				if($posStr && $opt{v}) {
					print "\t$posStr\tExon: ".get_string($exon)."\n";
				}
			}
		}
	}
}

my $fPtr;
if($fileOut) {
	if(!open($fPtr,">$fileOut")) {
		print "WARN: Could not open $fileOut for writing, using STDOUT\n";
	}
}
if(!$fPtr) {
	print "##SUMMARY\n\n";
}
foreach my $geneID (keys %targets) {
		my $chr = $targets{$geneID}{'chr'};
	foreach my $pos (@{$targets{$geneID}{'pos'}}) {
		my $ref = \%{$summary{$geneID}{$pos}};
		my $str;
#		$str .= "*" if $ref->{'*'};
#		$str .= "<" if $ref->{'<'};
#		$str .= ">" if $ref->{'>'};
#		$str .= "5" if $ref->{'5'};
#		$str .= "3" if $ref->{'3'};
		$str .= "5" if $ref->{'5'};
		$str .= "<" if $ref->{'<'};
		$str .= "*" if $ref->{'*'};
		$str .= ">" if $ref->{'>'};
		$str .= "3" if $ref->{'3'};
		if($fPtr) {
			print $fPtr "$geneID\t$chr\t$pos\t$str\n";
		} else {
			print "$geneID\t$chr\t$pos\t$str\n";
		}
	}
}
if($fPtr) {
	close($fPtr);
}
exit;

sub usage() {
	print "Syntax:\n$0\n";
	print << "EOF";
	-h      this message

	Either -l or -g need to be given
	[-l <file>]	tab-separated file with listing of queries.
			Either two formats:
				gene_name  chr   pos
				chr	   pos
	[-g <gene>]	just analyze this gene
	[-p <pos>]   	if a gene is given with -g, verify this coordinate
			It is only base position in chromosone, the chromosome
			is defined by the gene itself
	[-o <file>]	Output file, if not given, STDOUT is used along with messages
	[-e <# base>]	Internal threshold. Maximum distance between junction site
			and center of exon to call "<" or ">".
			(Think of -e as in -exon)
	[-i <# base>]	External threshold. Maximum distance outside the exon to
			call "5" or "3".
			(Think if -i as in -intron)
	[-v]		Verbose mode
EOF
exit;
}

sub mymax() {
	my @data = @{shift @_};
	my $max = 0;
	foreach my $item (@data) {
		$max = $item unless $item < $max;
	}
	return $max;
}

sub get_string {
	my $feature = shift;
	my $short;
	if(scalar(@_)) {
		$short = shift;
	}
	my $stable_id = $feature->stable_id;
	if($short) {
			$stable_id="";
	}
	my $seq_region = $feature->slice->seq_region_name;
	my $start = $feature->start;
	my $end = $feature->end;
	my $strand = $feature->strand;
	return "$stable_id $seq_region:$start-$end($strand)";
}

