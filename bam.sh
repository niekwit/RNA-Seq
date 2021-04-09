#!/bin/bash

threads=40
index_path="/home/niek/Documents/references/hisat2-index/GRCh37-hg19ucsc/hg19-index"


if [[ ! -d ./trim_galore ]]; then
	for file in fastq/*_1.fastq.gz 
	do
		echo "Trimming fastq files"
		file2="${file%_1.fastq.gz}_2.fastq.gz"
		trim_galore -j 4 -o ./trim_galore --paired $file $file2 2>> trim.log
	done
else
		echo "Skipping trimming (already performed)"
fi

echo "Aligning fastq files"
	mkdir -p bam
	for file in trim_galore/*1_val_1.fq.gz
	do
		echo $file
		hisat2_output="${file%1_val_1.fq.gz}-sort-bl.bam"
		hisat2_output="bam/${hisat2_output##*/}"
		file2="${file%1_val_1.fq.gz}2_val_2.fq.gz"
		hisat2 -p "$threads" -x "$index_path" -1 $file -2 $file2 2>> align.log | samtools view -q 15 -F 260 -bS -@ "$threads" - | samtools sort -@ "$threads" - > "$hisat2_output" #bam file will not contain unmapped and multimapped reads 
		samtools index -@ $threads $hisat2_output
	done

bam_list=$(ls bam/*.bam)

#create compressed numpy array for PCA
multiBamSummary bins -p $threads --bamfiles $bam_list -o multiBamSummary.npz

#PCA
plotPCA -in multiBamSummary.npz --colors #000000 #000000 #A9A9A9 #A9A9A9 #FF1493 #FF1493 #DDA0DD #DDA0DD #00008B #00008B #40E0D0 #40E0D0 #228B22 #228B22 #FF0000 #FF0000 #FFA500 #FFA500 --outFileNameData PCAdata.tab -o PCAplot.pdf
