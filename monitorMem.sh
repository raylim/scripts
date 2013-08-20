#!/bin/bash
# Serial version writes snapshots to $OUT.$i

OUT=$1

if [ ! $OUT ]; then
    OUT=/dev/stdout
fi

i=0
for i in `seq 1 40`; do
    ip=`printf "%02d" $i`
    free -m > $OUT.$ip
    echo "==========================================================" >> $OUT.$ip
    top -b -n 2 >> $OUT.$ip
    sleep 1
done
