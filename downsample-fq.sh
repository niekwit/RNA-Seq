#!/bin/bash

#WORKIND_DIR=$(pwd)
sample_size=$1

for read1 in *R1.fastq.gz
do
	SEED=$RANDOM
	read1_out="downsampling/${read1##*/}"
	read2="${read1%R1.fastq.gz}R2.fastq.gz"
	read2_out="downsampling/${read2##*/}"
	echo "Downsampling $read1 to $sample_size reads"
	seqtk sample -s $SEED  <(zcat $read1) $sample_size | gzip > $read1_out
	echo "Downsampling $read2 to $sample_size reads"
	seqtk sample -s $SEED  <(zcat $read1) $sample_size | gzip > $read2_out
done
