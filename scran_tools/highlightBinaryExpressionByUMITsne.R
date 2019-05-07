args <- commandArgs()

help <- function(){
  cat("highlightBinaryExpressionByUMITsne.R :
      - Highlight genes in tsne plots by whether they are expressed or not.
      - Colors must be brewer pal\n")
  cat("Usage: \n")
  cat("--RDS     : rds list file with phenodata, counts, and featureData               [ required ]\n")
  cat("--outDir : Output Directory                                                     [ required ]\n")
  cat("--Genes   : Common gene names to highlight in tSNE                              [ required ]
      (must be separated by a , and in the featureData table)               \n")
  cat("--Thresh  : Number of UMIs above which a gene is considered to be 'expressed'   [ default = 1 ]\n")
  cat("--Cols   : brewer pal palette                                                   [ default = Reds ]\n")
  cat("\n")
  q()
}

## Save values of each argument
if(!is.na(charmatch("--help",args)) || !is.na(charmatch("-h",args)) ){
  help()
} else {
  RDS     <-sub('--RDS=', '', args[grep('--RDS=', args)])
  outDir  <- sub( '--outDir=', '',args[grep('--outDir=',args)])
  Genes   <- sub( '--Genes=', '',args[grep('--Genes=',args)])
  Thresh  <- sub( '--Thresh=', '',args[grep('--Thresh=',args)])
  Cols    <- sub( '--Cols=', '',args[grep('--Cols=',args)])
}

# defaults
outDir <- paste(outDir, "/", sep="")
dir.create(outDir,recursive=TRUE, showWarnings = FALSE)

if (identical(Cols,character(0))){
  Cols <- "Reds"
}

if (identical(Thresh,character(0))){
  Thresh <- 1
}

library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(cowplot)

RDS <- "../non-normalized/P4136/data/Robjects/ExpressionList_QC_norm_clustered.rds"
outDir <- "Test"
Genes <- "SOX2,MYC"
Cols <- "Reds"
Thresh <- 1

# Load Data
rnd_seed <- 300
dataList <- readRDS(RDS)
m        <- dataList$counts
pD       <- dataList$phenoData
fD       <- dataList$featureData

m  <- m[,pD$PassAll]
pD <- pD[pD$PassAll,]

# set genes and titles
titles <- genes <- strsplit(Genes, split='\\,')[[1]]
genes <- gsub("\n","",genes)

# Subset m to look at raw UMI counts of GOI as a parameter for expression of that gene
m.sub <- t(m)
colnames(m.sub) <- sub("-", ".", fD$symbol)
genes <- genes[genes %in% colnames(m.sub)]

m.sub <- m.sub[,genes]
m.sub <- as.data.frame(m.sub)
m.sub[m.sub >= Thresh] <- 1

if (length(genes)==1) {
  names(m.sub) = genes[1]
}

m.sub$barcode <- rownames(m.sub)
stopifnot(m.sub$barcode==pD$barcode)
forPlot <- merge(pD, m.sub)
pal     <- colorRampPalette(brewer.pal(n=7,name=Cols))( length(unique(forPlot$Condition)) )

head(forPlot)
plots <- list()

for (i in seq_along(genes)) {
  gene <- genes[i]
  print(gene)
  p <- ggplot(forPlot, aes_string(x="tSNE1", y="tSNE2", colour=gene)) +
    geom_point(size=1) +
    scale_color_gradientn(colors=pal) +
    ggtitle(paste("Aggregated tSNE for", titles[i])) +
    theme_void(base_size=08) +
    theme(plot.title=element_text(size=rel(1)))
  plots[[gene]] <- p
}

plots_clust <- list()

for (i in seq_along(genes)) {
  gene <- genes[i]
  print(gene)
  p <- ggplot(forPlot, aes_string(x="tSNE1", y="tSNE2", colour=gene)) +
    geom_point(size=1) +
    facet_wrap(~Cluster) +
    scale_color_gradientn(colors=pal) +
    ggtitle("Faceted by cluster") +
    theme_void(base_size=08) +
    theme(plot.title=element_text(size=rel(1)))
  plots_clust[[gene]] <- p
}

plots_cond <- list()

for (i in seq_along(genes)) {
  gene <- genes[i]
  print(gene)
  p <- ggplot(forPlot, aes_string(x="tSNE1", y="tSNE2", colour=gene)) +
    geom_point(size=1) +
    facet_wrap(~Condition) +
    scale_color_gradientn(colors=pal) +
    ggtitle("Faceted by Condition") +
    theme_void(base_size=08) +
    theme(plot.title=element_text(size=rel(1)))
  plots_cond[[gene]] <- p
}

for (i in seq_along(genes)){
  gene <- genes[i]
  print(gene)
  pdf(sub("$", paste0("tSNE.", i,".UMI_Coloured.pdf"), outDir), 4, 4)
  print({
    plots[gene]
  })
  dev.off()
}

for (i in seq_along(genes)){
  gene <- genes[i]
  print(gene)
  pdf(sub("$", paste0("faceted.tSNE.", gene,".UMI_Coloured.pdf"), outDir), 4, 4)
  print({
    plot_grid(plots[[gene]],plots_clust[[gene]],plots_cond[[gene]])
  })
  dev.off()
}

