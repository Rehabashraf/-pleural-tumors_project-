#  Pleural Mesothelioma Transcriptomic Analysis using Limma-voom


## Project Overview 


This project performs a comprehensive bulk RNA-seq analysis of pleural mesothelioma histological subtypes using the Limma-voom workflow

 
The analysis includes quality control, normalization, differential expression analysis, pathway enrichment and Gene Set Enrichment Analysis (GSEA) 


to characterize transcriptomic differences among mesothelioma subtypes


## Objectives


. Perform differential expression analysis using Limma-voom.

. Compare transcriptomic profiles between pleural mesothelioma subtypes.

. Identify significantly differentially expressed genes.

. Explore enriched biological pathways using:

. Reactome

. GSEA


## Dataset

GEO accession : GSE327166


## Visualizations


The project generates:

 
. MDS Plot

. Density Plot

. Boxplots

. Sample Correlation Heatmap

. Hierarchical Clustering

. MA Plots

. Volcano Plots

. Heatmaps

. Venn Diagram

. GSEA Dotplots

. Ridgeplots

. Enrichment Maps


### Functional Analysis


Differential Expression


Genes were considered significant using:


Adjusted P-value < 0.05

|log2 Fold Change| > 1


## Gene Set Enrichment Analysis (GSEA)


GSEA was performed using ranked moderated t-statistics from Limma



### Reactome Pathway Analysis

Over-representation analysis (ORA) was performed using significant differentially expressed genes


## R Packages

. edgeR

. limma

. clusterProfiler

. ReactomePA

. enrichplot

. org.Hs.eg.db

. EnhancedVolcano

. Glimma

. ggplot2

. pheatmap

.dplyr


### Output Files


### Differential Expression Tables

epi_vs_mix.csv

epi_vs_Sar.csv

epi_vs_Mes.csv

mix_vs_Sar.csv

Sar_vs_Mes.csv

### Significant DEGs


sig_epi_vs_Sar.csv

sig_epi_vs_Mes.csv

sig_mix_vs_Sar.csv



### GSEA

gsea_epi_vs_mix.csv

gsea_epi_vs_Sar.csv

gsea_epi_vs_Mes.csv

gsea_mix_vs_Sar.csv

gsea_Sar_vs_Mes.csv

### Reactome

reactome_epi_vs_Sar.csv


