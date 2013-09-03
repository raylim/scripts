#!/bin/bash
# rename fastq files and create samples.txt and possibly samples.split.txt

prename 's/([^_]+)_[^_]+_L[^_]+_R([12])_([0-9]+)/$1_$3.$2/' fastq/*.fastq.gz
'ls' fastq/*.fastq.gz | sed 's:.*/::; s/_.*//' | sort | uniq > samples.txt
paste <('ls' fastq/*.fastq.gz | sed 's:.*/::; s/[._].*//') <('ls' fastq/*.fastq.gz | sed 's:.*/::; s/\..*//') | awk 'BEGIN { OFS = "\t" } $1 != $2 { print }' > samples.split.txt
if [ -z `cut -f1 samples.split.txt | uniq -d` ]; then
    prename 's/_.*//' fastq/*.fastq.gz
    rm samples.split.txt
fi
for x in `cut -f1 samples.split.txt | uniq -u`; do
    prename 's/_.*//' fastq/${x}_*.fastq.gz
    sed -i "/^$x\t/d" samples.split.txt
done
