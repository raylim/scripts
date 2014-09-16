cd /ifs/e63data/reis-filho/
find data projects -type d \
    \( -name bam -o -name tables -o -name alltables -o -name vcf \) \
    ! -path "*/log/*" ! -path "*/tmap/*" ! -path "*/gatk/*" ! -path "*/hydra/*" ! -path "*/bwa/*" \
    ! -path "*/varscan/*" ! -path "*/mutect/*" ! -path "*/scalpel/*" ! -path "*/som_sniper/*" -print0 | \
    rsync --verbose --progress -a -0 --files-from=- --prune-empty-dirs ./ /mount/limr/zedshared
