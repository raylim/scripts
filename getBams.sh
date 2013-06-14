#!/bin/bash
# usage: sh getBams.sh [-c] [-b] samples.txt

# $Id: $

SAMPLE_DIRS=$HOME/projects/sample_dirs.txt

copylink="ln -s"

# copy instead of link
if [ $1 == "-c" ]; then
    copylink="cp";
    shift;
fi

# re-build index
if [ $1 == "-b" ]; then
    find -L /projects/analysis /projects/seq_analysis /archive/analysis* /archive/solexa* -maxdepth 2 -name 'HS*' -or -name 'A*' > $SAMPLE_DIRS;
    shift;
fi


mkdir -p gsc_bam;
for sample in `cat $1`; do
    dir=`grep "$sample" $SAMPLE_DIRS`;
    bams=`find $dir -name "*dupsFlagged.bam"`
    numbams=`echo $bams | wc -w`
    if [ $numbams -eq 0 ]; then
        #echo "Found multiple possible bams for $sample: $bam";
        bams=`find $dir -name "*.bam"`
        numbams=`echo $bams | wc -w`;
    fi
    if [ $numbams -gt 1 ]; then
        #echo "Found multiple possible bams for $sample: $bam";
        bams=`echo $bams | grep "meta_bwa"`;
        numbams=`echo $bams | wc -w`;
    fi
    if [ $numbams -gt 1 ]; then
        #echo "Found multiple possible bams for $sample: $bam";
        max_lanes=`echo $bams | sed 's/ /\n/' | sed 's/.*\([0-9]\+\)_lanes.*/\1/' | awk 'BEGIN {max = 0} {if ($1 > max) max = $1} END {print max}'`;
        bams=`echo $bams | grep "${max_lanes}_lanes"`;
        numbams=`echo $bams | wc -w`;
    fi
    if [ $numbams -gt 1 ]; then
        echo "WARNING: Found > 1 bam for $sample in $dir";
    fi
    if [ $numbams -ge 1 ]; then
        bam=`echo $bams | cut -f 1 -d" "`
        #echo "Found $sample in $dir";
        $copylink -v $bam gsc_bam/$sample.bam;
        if [ -e $bam.bai ]; then
            $copylink $bam.bai gsc_bam/$sample.bam.bai;
        fi
    else
        echo "ERROR: Unable to find $sample in $dir";
    fi
done
