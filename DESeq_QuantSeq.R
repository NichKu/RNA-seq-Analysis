###################################################################################################
#
#     Set up working environment and import data
#
###################################################################################################

library(tximport)
library(RColorBrewer)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(vsn)
library(biomaRt)
library(DEGreport)

#library(readr)
#library(EnsDb.Hsapiens.v86)
#library(pathview)
#library(reshape2)
#library(tidyverse)
#library(gplots)
#library(DOSE)
#library(pathview)
#library(AnnotationHub)
#library(ensembldb)
#library(ggnewscale)


## Change to location of results directory with the amples.txt file and sample directories containing the quant.sf files 
wddir = "/Users/nicholas/Desktop/Results_quant_sf_DESeq/"
setwd(wddir)
base_dir = getwd()

# Create directory for output results
dir.create(file.path(getwd(), 'DESeq_output'), showWarnings = FALSE)

## Import sample and condition file
samples = read.table(file.path(base_dir, "samples.txt"), header = TRUE, stringsAsFactors=FALSE)
samples$condition <- factor(samples$condition)
samples$patient <- factor(samples$patient)
#samples$patient <- relevel()
samples        # Prints the sample / condition list


## Import Salmon quant files and create counts table
files <- file.path(base_dir, samples$sample, "quant.sf")
files
all(file.exists(files))        # Verify names of files match names in samples.csv, should return True
names(files)=samples$sample
samples$sample
names(files)
head(files)

rowdata = read.csv(file.path(base_dir, "salmon_tx2gene.tsv"), sep="\t", header = FALSE)
colnames(rowdata) = c("tx", "gene_id", "gene_name")
tx2gene = rowdata[,1:2]
tx2gene

txi <- tximport(files, 
                type = "salmon", 
                tx2gene = tx2gene)


head(txi$counts)               # This is the counts table for all of our samples

## Now to import the data into a DESeq Data Set (dds)
## Verify that sample names and colnames are the same
identical(samples$sample,colnames(txi$counts))

## Create a DEseqDataSet from txi count table
dds <- DESeqDataSetFromTximport(txi, samples, ~ condition)
#dds <- estimateSizeFactors(dds)
rownames(samples) <- colnames(txi$counts)

resultsNames(dds)


###################################################################################################
#
#    EXPLORATORY DATA ANALYSIS
#
###################################################################################################

## Set color palette for figures
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

## Perform a rlog transformation on count data (essentially a puts on a log2 scale)
## This helps our data assume a normal distribution and is good to do before these analyses
rld <- rlog(dds, blind=TRUE) # apply a regularized log transformation, ignoring information about experimental groups

## Setup annotation file to show the conditions on the figures
samples
treat_ann <- samples[,c("condition", "patient")]
treat_ann
#rownames(treat_ann) <- treat_ann$sample
#treat_ann$sample <- NULL
#treat_ann


## SAMPLE TO SAMPLE DISTANCE & CORRELATION HEATMAPS

## Sample correlation heatmap
corr_samps <- cor(as.matrix(assay(rld)))      # Computes pairwise correlations between samples based on gene expression
corr_samps
png(filename="DESeq_output/DESeq_sampleCorr_HM_sal.png", units = 'in', width = 12, height = 8, res = 250)
pheatmap(corr_samps,
         annotation = treat_ann,
         col=colors,
         main="Sample Correlations")
dev.off()

# Sample distance heatmap
sampleDists <- dist(t(assay(rld)))            # Computes Euclidean distance between samples based on gene expression
sampleDistMatrix <- as.matrix(sampleDists)

png(filename="DESeq_output/DESeq_sampleDist_HM_sal.png", units = 'in', width = 12, height = 8, res = 250)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         annotation = treat_ann,
         col=colors,
         main="Sample to Sample Distances")

dev.off()

## Principal Component Analysis
## Separates samples based on variation between sample's gene expression
## Greater variation will affect separation to a greater degree

data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)

percentVar <- round(100 * attr(data, "percentVar"))

png('DESeq_output/DESeq_PCA_sal.png', units='in', width=8, height=6, res=250)
ggplot(data, aes(PC1, PC2, color=condition)) +
  geom_point(size=3.5) +
  geom_text_repel(aes(label=name)) +
  scale_colour_manual(values = c("orange", "steelblue", 'red')) +
  theme_bw() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("")
dev.off()

###################################################################################################
#
#     DIFFERENTIAL EXPRESSSION ANALYSIS 
#
###################################################################################################
## DESeq = fx to calculate DE
## Combines multiple steps from DESeq
dds <- DESeq(dds, betaPrior = TRUE)


resultsNames(dds)
## QC
png(filename="DESeq_output/DispEst_sal.png", units = 'in', width = 12, height = 8, res = 250)
plotDispEsts(dds)
dev.off()

ntd <- normTransform(dds)
png(filename="DESeq_output/SDplot_sal.png", units = 'in', width = 12, height = 8, res = 250)
meanSdPlot(assay(ntd), ranks=F)
dev.off()

x <- assay(ntd)[,1]
y <- assay(ntd)[,2]
plot(.5*(x + y), y - x,
     cex=.5, col=rgb(0,0,0,.1), pch=20)
abline(h=0, col="red", lwd=3)

## Create full report
degResults(dds = dds, name = "pre versus post", org = NULL,
                        do_go = FALSE, group = "condition", xs = "condition",
                        path_results = "DESeq_output/")
resreport
resultsNames(dds)

##Highly Expressed Genes
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:100]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
png(filename="DESeq_output/Heatmap_highEG_sal.png", units = 'in', width = 12, height = 8, res = 250)
heatmap.2(assay(rld)[select,], col = hmcol, trace="none", margin=c(10, 6),
          labCol=colnames(dds), cexRow = 0.4)
dev.off()

## DE analysis
##uncomment if filtering is desired
#keep <- rowSums(counts(dds)) >= 0
#dds <- dds[keep,]
res <- results(dds, contrast = c("condition","post","pre"))
summary(res)
write.table(counts(dds),"DESeq_output/DESeq2.counts.tsv", sep="\t", quote=FALSE, col.names=NA)
write.table(counts(dds, normalized=T),"DESeq_output/DESeq2.counts_normalized.tsv", sep="\t", quote=FALSE, col.names=NA)

## Summary
hist (res$pvalue, breas=20)
hist (res$padj, breas=20)

## plot log2 fold changes
png(filename="DESeq_output/Log2_fold_sal.png", units = 'in', width = 12, height = 8, res = 250)
plotMA(res, ylim=c(-3,3))
dev.off()


## plot shrunken log2 fold changes using apeglm algorithm
#cannot be performed with betaPrior==True
#resLFC <- lfcShrink(dds, coef="condition_pre_vs_post", type="apeglm")
#summary(resLFC)

#png(filename="DESeq_output/shrunkLog2_fold.png", units = 'in', width = 12, height = 8, res = 250)
#plotMA(resLFC, ylim=c(-3,3))
#dev.off()


#resOrdered <- res[order(res$pvalue),]
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


## Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

## extract significant genes
threshold <- res$padj < padj.cutoff & abs(res$log2FoldChange) > lfc.cutoff
res$threshold <- threshold 
sigOE <- data.frame(subset(res, threshold==TRUE))
sigOE

sigOE_ordered <- sigOE[order(sigOE$padj), ]
top20_sigOE_genes <- rownames(sigOE_ordered[1:30, ])
top20_sigOE_genes



### use melt to modify the format of the data frame
#melted_top20_sigOE <- data.frame(melt(top20_sigOE_norm))

## check the column header in the "melted" data frame
#View(melted_top20_sigOE)

### add column names that make sense
#colnames(melted_top20_sigOE) <- c("gene", "samplename", "normalized_counts")


# Volcano plot
resOE_df <- data.frame(res)
resOE_df_ordered <- resOE_df[order(resOE_df$padj), ] 
resOE_df_ordered$genelabels <- rownames(resOE_df_ordered) %in% rownames(resOE_df_ordered[1:20,])

png(filename="DESeq_output/VolcanoPlot_annot_sal.png", units = 'in', width = 12, height = 8, res = 250)
ggplot(resOE_df_ordered) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = ifelse(genelabels == T, rownames(resOE_df_ordered),""))) +
  ggtitle("") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
dev.off()

### Annotate our heatmap (optional)
#annotation <- data.frame(sampletype=meta[,'sampletype'], 
#                         row.names=rownames(meta))
normalized_counts <- counts(dds, normalized=T)
norm_OEsig <- normalized_counts[rownames(sigOE),]
### Set a color palette
heat.colors <- brewer.pal(6, "YlOrRd")

norm_OEsig


## Setup annotation file to show the conditions on the figure
treat_ann_2 <- treat_ann[,c("condition", "patient")]
treat_ann_2

Var1        <- c("navy", "darkgreen")
names(Var1) <- c("post", "pre")
anno_colors <- list(Var1 = Var1)

### Run pheatmap
png(filename="DESeq_output/Heatmap_significantDEGs_sal.png", units = 'in', width = 12, height = 8, res = 250)
pheatmap(norm_OEsig, 
         color = heat.colors, 
         cluster_rows = T,
         show_rownames=T,
         border_color=NA, 
         fontsize = 10,
         scale="row",
         fontsize_row = 10, 
         height=20)
summary(res)


dev.off()
dev.new()
###################################################################################################
#
#     GO Enrichment Analysis
#
###################################################################################################

# Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genes <- as.character(rownames(res))
# Extract significant results
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ] 
signif_genes <- as.character(rownames(signif_res))
signif_genes
#keytypes(org.Hs.eg.db)
# Run GO enrichment analysis
ego <- enrichGO(gene = signif_genes,
                universe = all_genes,
                keyType = "SYMBOL",
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)


# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
png(filename="DESeq_output/dotplot_enrichGO.png", units = 'in', width = 12, height = 8, res = 250)
dotplot(ego, showCategory=50)
dev.off()


x2 <- pairwise_termsim(ego)
png(filename="DESeq_output/emaplot_enrichGO.png", units = 'in', width = 12, height = 8, res = 250)
emapplot(x2, showCategory=50)
dev.off()

# To color genes by log2 fold changes 
signif_res_lFC <- signif_res$log2FoldChange
png(filename="DESeq_output/cnetplot_enrichGO.png", units = 'in', width = 12, height = 8, res = 250)
cnetplot(ego,
         categorySize="pvalue",
         showCategory = 5,
         foldChange= signif_res_lFC,
         vertex.label.font=6)


dev.off()
###################################################################################################
#
#     Gene Set Enrichment Analysis 
#
###################################################################################################

mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(filters="hgnc_symbol", attributes=c("hgnc_symbol", "entrezgene_id"), values= all_genes, mart=mart)

indNA = which(is.na(genes$entrezgene_id))
genes_noNA <- genes[-indNA,]
indnodup = which(duplicated(genes_noNA$entrezgene_id) == F) 
genes_noNA_nodup <- genes_noNA[indnodup,]

lFC <- res$log2FoldChange[-indNA]
lFC <- lFC[indnodup]

names(lFC) <- genes_noNA_nodup$entrezgene_id
lFC

# Sort fold changes in decreasing order
lFC <- sort(lFC, decreasing = TRUE)
lFC

gseaKEGG <- gseKEGG(geneList = lFC,
                    organism = "hsa",
                    nPerm = 1000, # default number permutations
                    minGSSize = 5, # minimum gene set size
                    pvalueCutoff = 0.05, # padj cutoff value
                    verbose = FALSE)

# Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result
gseaKEGG_results




#Compared RNAseq to Nanostring
