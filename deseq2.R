library(tximport)
library(readr)
library(DESeq2)

args <- commandArgs(trailingOnly=TRUE)
dir <- args[1]
gtf <- args[2]

gtf <- "/home/niek/Documents/references/gtf/Homo_sapiens.GRCh38.95.gtf"
dir <- "/home/niek/Documents/analyses/RNA-Seq/40-449672264/delivery/fastq/salmon"

#Create sample list
samples <- read.table(file.path(dir, "samples.txt"), header=TRUE)
files <- file.path(dir,samples$dir, "quant.sf")
names(files) <- paste0("sample", 1:length(files))

#Create txdb from GTF file if it does not exist
txdb.filename <- paste0(gsub(".gtf","",gtf, fixed=TRUE),".txdb")

if(file.exists(txdb.filename) == FALSE){
  library(GenomicFeatures)
  txdb <- makeTxDbFromGFF(gtf)
  saveDb(txdb, txdb.filename)
  txdb <- loadDb(txdb.filename)
}

#Import transcript-level estimates
k <- keys(txdb,keytype="TXNAME")
tx2gene <- select(txdb,k,"GENEID","TXNAME")
txi <- tximport(files,type="salmon",tx2gene=tx2gene,ignoreTxVersion=TRUE)

