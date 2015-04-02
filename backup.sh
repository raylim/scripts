#!/bin/bash

PID=~/.backup.pid
TMP=`mktemp`;
if mountpoint -q "/mount/limr/zedshared/"; then
    if [[ -s $PID ]] || [[ -z `pgrep -F $PID | xargs ps` ]]; then 
        echo "searching for files in /ifs/e63data/reis-filho/"
        cd /ifs/e63data/reis-filho/ && \
        find data projects -type d \
            \( -name bam -o -name tables -o -name alltables -o -name vcf \) \
            ! -path "*/log/*" ! -path "*/tmap/*" ! -path "*/gatk/*" ! -path "*/hydra/*" ! -path "*/bwa/*" \
            ! -path "*/varscan/*" ! -path "*/mutect/*" ! -path "*/scalpel/*" ! -path "*/som_sniper/*" ! -path "*/rawdata/*" \
            ! -path "*/unprocessed_bam/*" ! -path "*/defuse/*" ! -path "*/chimscan/*" -print0 > ${TMP} &
        echo $! > $PID
        wait

        while [ 1 ]; do
            cd /ifs/e63data/reis-filho/ && \
            rsync --verbose --checksum --progress --stats --recursive --append --partial -a -0 --files-from=${TMP} --prune-empty-dirs ./ /mount/limr/zedshared &
            echo $! > $PID
            wait
            if [ "$?" = "0" ]; then
                echo "rsync complete"
                exit
            else
                echo "rsync failure, retrying in 1 minute..."
                sleep 60
            fi
        done
        rm $PID
        rm ${TMP}
    fi
fi
