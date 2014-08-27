#!/bin/bash
# rename fastq files and create samples.txt and possibly samples.split.txt

for x in fastq/*.fastq.gz; do
    if [ `grep -o "_" <<< "$x" | wc -l` -gt 1 ]; then
        prename 's/-//g; s/(.+)_([ATGC]{6,8})_L([^_])+_R([12])_([0-9]+)/$1-$2$3$5.$4/; s/(.+)_[^_]+_L([^_])+_R([12])_([0-9]+)/$1_$2$4.$3/; s/_//g; s/-/_/g' $x;
    fi;
done
'ls' fastq/*.fastq.gz | sed 's:.*/::; s/_.*//' | sort | uniq > samples.txt
paste <('ls' fastq/*.fastq.gz | sed 's:.*/::; s/[._].*//') <('ls' fastq/*.fastq.gz | sed 's:.*/::; s/\..*//') | awk 'BEGIN { OFS = "\t" } $1 != $2 { print }' | uniq > samples.split.txt
#if [ -z `cut -f1 samples.split.txt | uniq -d` ]; then
#    prename 's/_.*//' fastq/*.fastq.gz
#    rm samples.split.txt
#fi
#for x in `cut -f1 samples.split.txt | uniq -u`; do
#prename 's/_.*//' fastq/${x}_*.fastq.gz
#sed -i "/^$x\t/d" samples.split.txt
#done
