#!/bin/bash
# usage sh findBam.sh sample_dirs.txt {SAMPLE ID}
# $Id: $

INDEX_SAMPLE_DIRS=~/share/scripts/indexSampleDirs.sh

if [ $1 == "-b" ]; then
    echo "Indexing...";
    sh ${INDEX_SAMPLE_DIRS} $1;
    shift;
fi

sample_dirs=$1
sample=$2
dir=`grep "$sample" $sample_dirs`;
if [ "$dir" = "" ]; then
    echo "Unable to find sample: $sample" >&2;
    exit 1;
fi
bams=`find $dir -name "*dupsFlagged.bam"`
numbams=`echo $bams | wc -w`
if [ $numbams -eq 0 ]; then
    bams=`find $dir -name "*.bam"`
    numbams=`echo $bams | wc -w`;
fi
if [ $numbams -gt 1 ]; then
    echo "Found multiple possible bams for $sample: $bams" >&2;
    bams=`echo $bams | tr ' ' '\n' | grep "meta_bwa" | grep -v "old"`;
    numbams=`echo $bams | wc -w`;
fi
if [ $numbams -gt 1 ]; then
    echo "Found multiple possible bams for $sample: $bams" >&2;
    max_lanes=`echo $bams | sed 's/ /\n/' | sed 's/.*\([0-9]\+\)_lanes.*/\1/' | awk 'BEGIN {max = 0} {if ($1 > max) max = $1} END {print max}'`;
    bams=`echo $bams | tr ' ' '\n' |  grep "${max_lanes}_lanes"`;
    numbams=`echo $bams | wc -w`;
fi
if [ $numbams -gt 1 ]; then
    echo "WARNING: Found > 1 bam for $sample in $dir" >&2;
fi
if [ $numbams -ge 1 ]; then
    bam=`echo $bams | cut -f 1 -d" "`
    echo "Found $sample in $dir" >&2;
    echo $bam;
else
    echo "ERROR: Unable to find $sample in $dir" >&2;
    exit 1;
fi

