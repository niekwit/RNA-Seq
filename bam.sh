#!/bin/bash

threads=40
index_path="/home/niek/Documents/references/hisat2-index/GRCh37-hg19ucsc/hg19-index"

#transcript_index_path="/home/niek/Documents/references/fasta/Human/gencode.v35.transcripts/salmon_index"
#gtf="/home/niek/Documents/references/gtf/Homo_sapiens.GRCh38.95.gtf"

if [[ ! -d ./trim_galore ]]; then
	echo "Trimming fastq files"
	for read1 in *R1.fastq.gz 
	do
		
		echo $read1
		read2="${read1%_R1.fastq.gz}_R2.fastq.gz"
		trim_galore -j 4 -o ./trim_galore --paired $read1 $read2 2>> trim.log
	done
else
		echo "Skipping trimming (already performed)"
fi

echo "Aligning fastq files"
	mkdir -p bam
	for file in trim_galore/*1_val_1.fq.gz
	do
		echo $file
		echo $file >> align.log
		hisat2_output="${file%_R1_val_1.fq.gz}-sort-bl.bam"
		hisat2_output="bam/${hisat2_output##*/}"
		#salmon_output="salmon/${hisat2_output##*/}-quant"
		file2="${file%1_val_1.fq.gz}2_val_2.fq.gz"
		hisat2 -p "$threads" -x "$index_path" -1 $file -2 $file2 2>> align.log | samtools view -q 15 -F 260 -bS -@ "$threads" - | samtools sort -@ "$threads" - > "$hisat2_output" #bam file will not contain unmapped and multimapped reads 
		samtools index -@ $threads $hisat2_output
		
	done

bam_list=$(ls bam/*.bam)

#create compressed numpy array for PCA
multiBamSummary bins -p $threads --bamfiles $bam_list -o multiBamSummary.npz

#PCA
plotPCA -in multiBamSummary.npz --outFileNameData PCAdata.tab -T "PCA of BAM files" -o PCAplot.pdf --colors #ff0000,#ff0000,#0000ff,#0000ff,#3cb371,#3cb371,#ee82ee,#ee82ee,#ffa500,#ffa500,#8ae600,#8ae600,#404040,#404040,#a0a0a0,#a0a0a0,#7b3f26,#7b3f26


#ff0000,#0000ff,#3cb371,#ee82ee,#ffa500,#8ae600,#404040,#a0a0a0,#7b3f26
