#!/usr/bin/perl

$outdir=$ARGV[0];
my %snp=();

open(IN,"all_normal_SNPs.txt") || die "cannot open SNP file\n";
while(<IN>){
    chomp;
    my @fields = split(/\s+/,$_);
    my $chrpos = $fields[0];
    $snp{$chrpos} = $_;
}

@files=`ls $outdir/*_snvmix2_novel_codon_annot_artifact_filtered.txt`;

foreach $file (@files){
    chomp($file);
    ($name,$rest)=split(/\./,$file);
    $fileout=$name."_nfilter.txt";
    open(OUT,">$fileout") || die "cannot open $outdir/$fileout file\n";
    open(IN,"$file") || die "cannot open $file\n";
    while(<IN>){
	chomp;
	($chrpos,$ensid,$strand,$refcodon,$posmut,$refbase,$mutbase,$readsupport,$refcodon1,$mutcodon,$refaa,$mutaa,$change,$charge,$polarity,$score,$gene,$relposition)=split(/\s+/,$_);
	#print "$chrpos\n";	
	$a=$snp{$chrpos};
	if(!(defined($a))) {
	    print OUT "$_\n";
	}
	else{
	    #print "$chrpos $file\n";
	}
	
    }
}
	
