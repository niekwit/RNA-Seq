library(ggplot2)
library(ggrepel)

setwd("/home/niek/Documents/scripts/RNA-Seq")
df <- read.csv(file="PCAdata.csv")
df$replicate <- as.factor(df$replicate)

colours <- c("Brian-Control-Nx"="#ff0000",
             "Control-Hx"="#0000ff",
             "Control-Nx"="#3cb371",
             "HIF1B-Hx"="#ee82ee",
             "HIF1B-Nx"="#ffa500",
             "USP43-C19-Hx"="#8ae600",
             "USP43-C19-Nx"="#404040",
             "USP43mp-Hx"="#a0a0a0",
             "USP43mp-Nx"="#7b3f26")

p <- ggplot(df, aes(x = `PC1`, y = `PC2`)) +
  geom_point(aes(size = `Read.count`,shape = `replicate`, fill = `sample`),alpha=0.6) +
  theme_bw() +
  labs(size="mapped paired read counts") +
  scale_shape_manual(values=c(21,22)) +
  scale_fill_manual(values=c(colours)) + 
  scale_size(range = c(7, 17)) +
  xlab("PC1 (86.9% of var. explained)") +
  ylab("PC2 (8.3% of var. explained)") +
  guides(fill=guide_legend(override.aes=list(shape=21, size=5)),
         shape=guide_legend(override.aes=list(size=5)))

ggsave("PCAplot.pdf",plot=p, width=8,height=7) 

