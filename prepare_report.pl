# This perl script prepares a final SNVMix report 
# Sanja Rogic
#!/usr/bin/perl

open(OUT,">MCL_snv_final_report_filtered.txt");

@files=`ls *_snvmix2_novel_codon_annot_artifact_filtered_nfilter.txt`;

$print=join("\t","Sample","Gene Name","EnsembleID","Chromosome","Position","Strand","Reference codon","Position mutated","Mutated codon","Reference base","Mutated base","Read Support","Ref AA","Mutated AA","Change","Charge","Polarity","Score","Relative position");
print OUT "$print\n";
foreach $file (@files){
    chomp($file);
    ($name,$rest)=split("_",$file);
    open(IN,$file);
    while(<IN>){
	chomp;
	($chrpos,$ensid,$strand,$refcodon,$posmut,$refbase,$mutbase,$readsupport,$refcodon1,$mutcodon,$refaa,$mutaa,$change,$charge,$polarity,$score,$gene,$relposition)=split(/\s+/,$_);
	($chr,$pos)=split(":",$chrpos);
	$chr=~s/chr//;
	if($change eq "CODING"){
	    $print=join("\t",$name,$gene,$ensid,$chr,$pos,$strand,$refcodon,$posmut,$mutcodon,$refbase,$mutbase,$readsupport,$refaa,$mutaa,$change,$charge,$polarity,$score,$relposition);
	    print OUT "$print\n";
	}
	
    }
}
	
