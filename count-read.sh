#!/bin/bash

echo "Uniquely mapped read counts before deduplication:" >> mapped_read_count_no_dedup.txt	
for file in bam/*-sort-bl.bam
do
	count=$(samtools view -@ 40 -c -f 1 -F 12 $file)
	echo "$file: $count" >> mapped_read_count_no_dedup.txt
done
