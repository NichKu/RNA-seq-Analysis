**RNA-seq Data Analysis with DEGA, EA and GSEA**
---

author: Nicholas Küng\
contact: nicholas.kueng@extern.insel.ch


This script uses the Salmon quantification output, a samplesheet file and a gene annotation file to perform:
- Exploratory Data Analyis with a sample correlation and sample distance heatmap, and a PCA plot,
- Differential Expression Analysis,
- GO Enrichment Analysis and
- Gene Set Enrichment Analysis

To run the program and save the output define your working directory. The path has to be entered in line 34 as following:\
**wddir = "yourPATH/"**

In this directory the following files and subdirectories should exist:
- sample specific Salmon output files with quant.sf files
- sample.txt
- salmon_tx2gene.tsv

see [Salmon_QuantFiles](https://github.com/NichKu/RNA-seq-Analysis/tree/master/Salmon_QuantFiles)

**Dependencies**
- tximport
- RColorBrewer
- DESeq2
- ggplot2
- ggrepel
- pheatmap
- clusterProfiler
- enrichplot
- org.Hs.eg.db
- vsn
- biomaRt
- DEGreport