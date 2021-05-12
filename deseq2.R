suppressMessages(library(tximport))
suppressMessages(library(readr))
suppressMessages(library(DESeq2))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(biomaRt))
suppressMessages(library(ggplot2))


args <- commandArgs(trailingOnly=TRUE)
work.dir <- args[1]
gtf <- args[2]
script.dir <- args[3]
species <- args[4]

#Create sample files for all comparisons
samples.master <- read.table(file.path(work.dir,"samples.txt"), header=TRUE)
number.exp <- ncol(samples.master)-2
exp.names <- colnames(samples.raw[,3:ncol(samples.raw)])

df.list <- list()#list for storing comparison dfs
for (i in 1:length(exp.names)) { #create data frames for each experiment from master sample sheet
  exp <- exp.names[i]
  df.temp <- samples.master
  temp.list <- exp.names
 
  #select which columns to remove
  to.remove <- exp != temp.list
  drop.columns <- temp.list[to.remove]
  
  #create data frame with only one experiment (exp)
  df.temp <- df.temp[, !(names(df.temp) %in% drop.columns)]
  
  #remove NAs (keeps relevant experimental settings)
  df.temp <- df.temp[complete.cases(df.temp), ]
  
  #add df to df.list
  df.list[[i]] <- df.temp
  }

names(df.list) <- exp.names #name dfs in list

#Generate DESeq2 required variables
#Create txdb from GTF file if it does not exist
txdb.filename <- paste0(gsub(".gtf","",gtf, fixed=TRUE),".txdb")

if(file.exists(txdb.filename) == FALSE){
  txdb <- makeTxDbFromGFF(gtf)
  saveDb(txdb, txdb.filename)
  txdb <- loadDb(txdb.filename)
} else {txdb <- loadDb(txdb.filename)}

#Create transcript to gene file
k <- keys(txdb,keytype="TXNAME")
tx2gene <- select(txdb,k,"GENEID","TXNAME")

#select ensembl data base for ensemble ID to gene symbol conversion
if(species == "mouse"){
  library="EnsDb.Mmusculus.v79"
} else if(species == "human"){
  library="EnsDb.Hsapiens.v79"}
library(library)

#Run DESeq2 for each sample df in df.list
for (i in 1:length(df.list)){
  #create file list for DESeq2 input
  samples <- df.list[[i]]
  files <- file.path(work.dir,"salmon",paste0(samples$sample,"-quant"), "quant.sf")
  names(files) <- samples$sample
  
  txi <- tximport(files,
                  type="salmon",
                  tx2gene=tx2gene,
                  ignoreTxVersion=TRUE)
  
  #Construct DESeq2 data set
  dds <- DESeqDataSetFromTximport(txi,
                                  colData=samples,
                                  design= ~ condition)
  
  #Pre-filtering: remove rows with very few reads
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  #Set reference level: condition that is marked control in exp column
  ref <- which(samples == "control", arr.ind=TRUE) #select rows that contain control
  ref <- ref[,"row"]#get indeces of rows that contain control
  ref <- samples[ref[[1]],] #subset samples for row that contains just control condition
  ref <- ref[["condition"]] #extract reference condition
  dds$condition <- relevel(dds$condition, ref=ref)
  
  #Differential expression analysis
  dds <- DESeq(dds)
  res <- results(dds, alpha=0.01) #adjusted p-value < 0.01
  res <- res[order(res$padj),]
  
  #create directory for output
  dir.out <- paste0(work.dir,"/",resultsNames(dds)[2])
  if(dir.exists(dir.out) == FALSE){
    dir.create(dir.out)
  }
  
  #Generate MA plot
  resLFC <- lfcShrink(dds, 
                      coef=resultsNames(dds)[2], 
                      type="apeglm") #results for plotting (shrinkage of size effect)
  
  pdf(file=paste0(dir.out,"/MA-plot.pdf"))
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
  
  ggsave(paste0(dir.out,"/PCA-plot.pdf"),p)
  
  #Include normalised read counts in output table
  df <- as.data.frame(res)
  df$GENEID <- rownames(df)
  
  df.reads <- as.data.frame(counts(dds,normalized=TRUE)) #gets normalised read counts from dds
  df.reads$GENEID <- row.names(df.reads)
  df <- merge(x=df,
              y=df.reads,
              by="GENEID")
  
  #Convert Ensembl gene IDs to gene symbols
  genes <- res@rownames
  gene.symbols <- ensembldb::select(library, 
                                    keys= genes, 
                                    keytype = "GENEID", 
                                    columns = c("SYMBOL","GENEID"))
  
  suppressMessages(library(dplyr))
  df <- merge(x=df,
              y=gene.symbols,
              by="GENEID")
  
  #order data for padj
  df <- df[order(df$padj), ]
  
  write.csv(df, 
            file=paste0(dir.out,"/DESeq-output.csv"),
            row.names=FALSE)
  
}

