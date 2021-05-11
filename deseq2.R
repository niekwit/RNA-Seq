suppressMessages(library(tximport))
suppressMessages(library(readr))
suppressMessages(library(DESeq2))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(biomaRt))
library(ggplot2)
library(EnsDb.Hsapiens.v79)

args <- commandArgs(trailingOnly=TRUE)
work.dir <- args[1]
gtf <- args[2]
script.dir <- args[3]

#Create sample list
samples <- read.table(file.path(work.dir,"samples.txt"), header=TRUE)
files <- file.path(work.dir,"salmon",paste0(samples$sample,"-quant"), "quant.sf")
names(files) <- samples$sample

#Create txdb from GTF file if it does not exist
txdb.filename <- paste0(gsub(".gtf","",gtf, fixed=TRUE),".txdb")

if(file.exists(txdb.filename) == FALSE){
  txdb <- makeTxDbFromGFF(gtf)
  saveDb(txdb, txdb.filename)
  txdb <- loadDb(txdb.filename)
} else {txdb <- loadDb(txdb.filename)}

#Import transcript-level estimates
k <- keys(txdb,keytype="TXNAME")
tx2gene <- select(txdb,k,"GENEID","TXNAME")
txi <- tximport(files,type="salmon",tx2gene=tx2gene,ignoreTxVersion=TRUE)

#Construct DESeq2 data set
dds <- DESeqDataSetFromTximport(txi,
                                   colData=samples,
                                   design= ~ condition)

#Pre-filtering: remove rows with very few few reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Set reference level
dds$condition <- relevel(dds$condition, ref="normoxia")

#Differential expression analysis
dds <- DESeq(dds)
res <- results(dds, alpha=0.01) #adjusted p-value < 0.01
res <- res[order(res$padj),]

#Generate MA plot
resLFC <- lfcShrink(dds, 
                    coef=resultsNames(dds)[2], 
                    type="apeglm") #results for plotting (shrinkage of size effect)

pdf(file=paste0(work.dir,"/DESeq2/MA-plot.pdf"))
plotMA(resLFC, ylim=c(-2,2))
dev.off()

#Generate PCA plot
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition", "sample"), returnData=TRUE)
percentVar <- round(100* attr(pcaData, "percentVar"))
p <- ggplot(pcaData, aes(PC1, PC2, color=condition,shape=sample)) +
  theme_light() +
  theme(legend.position = "right") +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()

ggsave(paste0(work.dir,"/DESeq2/PCA-plot.pdf"),p)

#Convert Ensembl gene IDs to gene symbols
genes <- res@rownames
gene.symbols <- ensembldb::select(EnsDb.Hsapiens.v79, keys= genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
df <- as.data.frame(res)
df$GENEID <- rownames(df)

suppressMessages(library(dplyr))
df <- merge(x=df,y=gene.symbols,by="GENEID")

write.csv(df, 
          file=paste0(work.dir,"/DESeq2/DESeq-output.csv",rownames=FALSE))
