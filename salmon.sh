#!/bin/bash

threads=40
#index_path="/home/niek/Documents/references/hisat2-index/gencode.v35.transcripts/hisat2-index-transcripts-gencode-v35"

transcript_index_path="/home/niek/Documents/references/fasta/Human/gencode.v35.transcripts/salmon_index"
gtf="/home/niek/Documents/references/gtf/Homo_sapiens.GRCh38.95.gtf"

if [[ ! -d ./trim_galore ]]; then
	for read1 in *R1.fastq.gz
	do
		echo "Trimming fastq files"
		echo $read1
		read2="${read1%_R1.fastq.gz}_R2.fastq.gz"
		trim_galore -j 4 -o ./trim_galore --paired $read1 $read2 2>> trim.log
	done
else
		echo "Skipping trimming (already performed)"
fi

echo "Running Salmon"
	mkdir -p salmon
	for file in trim_galore/*1_val_1.fq.gz
	do
		echo $file
		echo $file >> salmon.log
		#hisat2_output="${file%_R1_val_1.fq.gz}-sort-bl.bam"
		salmon_output="salmon/${hisat2_output##*/}-quant"
		file2="${file%1_val_1.fq.gz}2_val_2.fq.gz"
		salmon quant --index $transcript_index_path -l A -g $gtf -p 12 -1 $file -2 $file2 --validateMappings --gcBias -o $salmon_output 2>> salmon.log
	done
