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

echo "Aligning fastq files"
	#mkdir -p bam
	mkdir -p salmon
	for file in trim_galore/*1_val_1.fq.gz
	do
		echo $file
		echo $file >> salmon.log
		hisat2_output="${file%_R1_val_1.fq.gz}-sort-bl.bam"
		#hisat2_output="bam/${hisat2_output##*/}"
		salmon_output="salmon/${hisat2_output##*/}-quant"
		file2="${file%1_val_1.fq.gz}2_val_2.fq.gz"
		#hisat2 -p "$threads" -x "$index_path" --rna-strandness FR -1 $file -2 $file2 2>> align.log | samtools view -q 15 -F 260 -bS -@ "$threads" - | samtools sort -@ "$threads" - > "$hisat2_output" #bam file will not contain unmapped and multimapped reads 
		#samtools index -@ $threads $hisat2_output
		salmon quant --index $transcript_index_path -l A -g $gtf -p 12 -1 $file -2 $file2 --validateMappings --gcBias -o $salmon_output 2>> salmon.log
	done

#bam_list=$(ls bam/*.bam)

#create compressed numpy array for PCA
#multiBamSummary bins -p $threads --bamfiles $bam_list -o multiBamSummary.npz

#PCA
#plotPCA -in multiBamSummary.npz --colors #000000 #000000 #A9A9A9 #A9A9A9 #FF1493 #FF1493 #DDA0DD #DDA0DD #00008B #00008B #40E0D0 #40E0D0 #228B22 #228B22 #FF0000 #FF0000 #FFA500 #FFA500 --outFileNameData PCAdata.tab -o PCAplot.pdf
