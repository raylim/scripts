#!/bin/bash

# $Id: $

FIND_BAM=~/share/scripts/findBam.sh
INDEX_SAMPLE_DIRS=~/share/scripts/indexSampleDirs.sh
SAMPLE_DIRS=$HOME/projects/sample_dirs.txt

if [ $1 == "-b" ]; then
    sh ${INDEX_SAMPLE_DIRS} ${SAMPLE_DIRS};
    shift;
fi

samples=`cat $1`;

mkdir -p bam;
cd bam;
for sample in ${samples}; do
    bam=`sh ${FIND_BAM} ${SAMPLE_DIRS} ${sample}`;
    ln -vs ${bam} ${sample}.bam;
done;

