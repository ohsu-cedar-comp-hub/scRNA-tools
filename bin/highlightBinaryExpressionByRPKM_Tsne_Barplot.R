args <- commandArgs()

help <- function(){
  cat("highlightBinaryExpressionByUMI_Tsne_BarPlot.R :
      - Highlight genes in tsne plots by whether they are expressed or not.
      - Colors must be brewer pal\n")
  cat("Usage: \n")
  cat("--RDS     : rds list file with phenodata, counts, and featureData               [ required ]\n")
  cat("--outDir  : Output Directory                                                    [ required ]\n")
  cat("--Genes   : Common gene names to highlight in tSNE                              [ required ]
      (must be separated by a , and in the featureData table)               \n")
  cat("--Thresh  : RPKM Threshold above which a gene is considered to be 'expressed'   [ default = 0.25 ]\n")
  cat("--Cols    : brewer pal palette                                                  [ default = Reds ]\n")
  cat("\n")
  q()
}

## Save values of each argument
if(!is.na(charmatch("--help",args)) || !is.na(charmatch("-h",args)) ){
  help()
} else {
  RDS     <- sub('--RDS=', '', args[grep('--RDS=', args)])
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
  Thresh <- 0.25
}

library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(cowplot)
library(edgeR)

# Genelength table
geneLength <- read.table(file="/home/groups/CEDAR/anno/CellRanger/GRCh38/refdata-cellranger-GRCh38-3.0.0.geneLength.txt.gz", 
                         stringsAsFactors = FALSE, header=TRUE)
head(geneLength)
rownames(geneLength) <- geneLength$id

# Load Data
rnd_seed <- 300
dataList <- readRDS(RDS)
m        <- dataList$counts
pD       <- dataList$phenoData
fD       <- dataList$featureData

m  <- m[,pD$PassAll]
pD <- pD[pD$PassAll,]

# Add gene length info to fD
iv <- match(rownames(m), rownames(geneLength))
head(geneLength[iv,])
head(fD)
rownames(fD) <- fD$id
stopifnot(rownames(geneLength[iv,])==rownames(fD))
fD$geneLength <- geneLength[iv,]$width
head(fD)

# set genes and titles
titles <- genes <- strsplit(Genes, split='\\,')[[1]]
genes <- gsub("\n","",genes)

# Generate RPKM matrix
m.rpkm <- rpkm(m, gene.length=fD$geneLength)
rownames(m.rpkm) <- sub("-", ".", fD$symbol)
genes <- genes[genes %in% rownames(m.rpkm)]

write.table(m.rpkm, file=paste(sub("$", "RPKM_counts_matrix.txt", Dir2)), quote = F, sep = "\t")

# log2 transform the normalized counts for genes of interest
expr <- t(m.rpkm)[,genes]

# combine the normalized expression with tSNE1 and tSNE2
add <- data.frame(expr,
                  barcode=colnames(m))

# Generate dataframe with all relevant expression information
## melt for ggplot
forPlot <- left_join(add,pD[,c("barcode","Condition")])
forPlot <- melt(forPlot,id=c("barcode","Condition")) %>%
  dplyr::rename(Expression=value) %>%
  dplyr::rename(Gene=variable)
## get box colors
pal     <- colorRampPalette(brewer.pal(n=7,name=Cols))( length(unique(forPlot$Condition)) )

forPlot$threshold_RPKM <- forPlot$Expression > Thresh

smryByGene <- group_by(forPlot, Gene) %>%
  summarise(NumCellsRPKM1 = sum(threshold_RPKM),
            TotalNumCells = n())
smryByGene$RPKMfraction = (smryByGene$NumCellsRPKM1) / (smryByGene$TotalNumCells)

smryByGene <- as.data.frame(smryByGene)
smryByGene <- merge(forPlot, as.data.frame(smryByGene))
smryByGene$Condition <- factor(smryByGene$Condition,levels=unique(smryByGene$Condition))

p1 <- ggplot(data=smryByGene, mapping=aes(x=Gene, y=RPKMfraction, fill=Condition)) +
  geom_bar(stat="identity") +
  ylab(paste("Number of Cells that meet RPKM >", Thresh)) +
  ggtitle("Cells expressing our GOI") +
  xlab("Gene") +
  scale_fill_manual(values=pal) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)
        ,panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))

pdf(sub("$","BarPlot_FractionCellsMeetRPKM_ByCondition.pdf",outDir), 4, 4)
print({
  p1
})
dev.off()

# Obtain information for colouring tSNE by presence/absence of gene expression
plots <- list()
forPlot$threshold_RPKM <- as.integer(as.logical(forPlot$threshold_RPKM))

# Split up by gene name and then run tSNE for each subsetted dataframe
for (i in seq_along(genes)) {
  gene <- genes[i]
  df <- filter(forPlot,Gene==gene)
  df <- df%>%select(-"barcode",everything())
  df <- left_join(df,pD[,c("barcode","tSNE1","tSNE2","Cluster")])
  p <- ggplot(df, aes(x=tSNE1, y=tSNE2, color=threshold_RPKM)) +
    geom_point(size=1) +
    scale_color_gradientn(colors=pal) +
    ggtitle(titles[i]) +
    theme_void(base_size=08) +
    theme(plot.title=element_text(size=rel(1)))
  plots[[gene]] <- p
}

plots_clust <- list()

for (i in seq_along(genes)) {
  gene <- genes[i]
  df <- filter(forPlot,Gene==gene)
  df <- df%>%select(-"barcode",everything())
  df <- left_join(df,pD[,c("barcode","tSNE1","tSNE2","Cluster")])
  p <- ggplot(df, aes(x=tSNE1, y=tSNE2, color=threshold_RPKM)) +
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
  df <- filter(forPlot,Gene==gene)
  df <- df%>%select(-"barcode",everything())
  df <- left_join(df,pD[,c("barcode","tSNE1","tSNE2","Cluster")])
  p <- ggplot(df, aes(x=tSNE1, y=tSNE2, color=threshold_RPKM)) +
    geom_point(size=1) +
    facet_wrap(~Condition) +
    scale_color_gradientn(colors=pal) +
    ggtitle("Faceted by Condition") +
    theme_void(base_size=08) +
    theme(plot.title=element_text(size=rel(1)))
  plots_cond[[gene]] <- p
}

for ( i in names(plots) ){
  print(i)
  pdf(sub("$", paste0("tSNE.", i, ".", Thresh, ".RPKM_Coloured.pdf"), outDir), 4, 4)
  print({
    plots[i]
  })
  dev.off()
}

for (i in names(plots)){
  print(i)
  pdf(sub("$", paste0("faceted.tSNE.", i,".RPKM_Coloured.pdf"), outDir), 4, 4)
  print({
    plot_grid(plots[[i]],plots_clust[[i]],plots_cond[[i]])
  })
  dev.off()
}

