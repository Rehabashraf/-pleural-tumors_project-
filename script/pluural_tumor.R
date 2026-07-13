#set working directory 

getwd()

setwd("/home/rehabashraf/R/x86_64-pc-linux-gnu-library/4.5/pleural_tumors/pleural_tumor2")

list.files()

# calling package
library(limma)
library(edgeR)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(RColorBrewer)
library(Glimma)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(EnhancedVolcano)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)

# load data

plura_count=read.delim("GSE327166_bulkRNA_rawcounts_matrix.tsv.gz",row.names = 1,header = TRUE,check.names = FALSE)

plura_count=as.matrix(plura_count)

storage.mode(plura_count)="integer"

plura_meta=read.delim("GSE327166_bulkRNA_sample_metadata.tsv.gz" ,header = TRUE,row.names = 1,check.names = FALSE)

# organise the data

colnames(plura_meta)[colnames(plura_meta)=="case_diagnosis"]="condition"

plura_meta$condition=as.character(plura_meta$condition)

plura_meta$condition[plura_meta$condition=="Mesothelioma-Epithelial"]="M_Epithelial"

plura_meta$condition[plura_meta$condition=="Mesothelioma-Mixed"]="M_mixed"

plura_meta$condition[plura_meta$condition=="Mesothelioma-Sarcomatoid"]="M_Sarcomatoid"

#subseting data

keep=plura_meta$condition %in% c("M_Epithelial","M_mixed","M_Sarcomatoid","Mesothelioma")

plura_meta=plura_meta[keep,]

plura_count=plura_count[,keep] 

# match data

all(rownames(plura_meta)==colnames(plura_count))

rownames(plura_meta)=colnames(plura_count)

all(rownames(plura_meta)==colnames(plura_count))

#convert category data to factor 

plura_meta$condition=as.factor(plura_meta$condition)

plura_meta$condition=factor(plura_meta$condition,levels = c("M_Epithelial", "M_mixed","M_Sarcomatoid","Mesothelioma"))

# creat an obeject

deg_object=DGEList(counts = plura_count,samples=plura_meta,group=plura_meta$condition)

# inspect DGEList

dim(deg_object)

deg_object$samples

# Add gene-level information 

deg_object$ENTREZ_ID=mapIds(org.Hs.eg.db,keys = row.names(plura_count),column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first")

deg_object$ENTREZ_ID

deg_object$symbol=mapIds(org.Hs.eg.db,keys = rownames(plura_count),column = "SYMBOL",keytype ="ENSEMBL",multiVals = "first" )

deg_object$symbol

#filter low experssion gene by filterbyexpxpr

filter_gene=filterByExpr(deg_object,group = plura_meta$condition, min.count=10,min.total.count=15,min.prop=0.7)

deg_object_filter=deg_object[filter_gene, ,keep.lib.siz=FALSE]

#print gete before and after filteration

cat("gene before filter:",nrow(deg_object),"\n")

cat("gene after filter:",nrow(deg_object_filter),"\n")

cat("generemove:",sum(!filter_gene),"\n")

table(filter_gene)

#TMM nornalization 

deg_norm=calcNormFactors(deg_object_filter,method = "TMMwsp" )

deg_norm$samples$norm.factors

#log cpm for visualization 

deg_norm_cpm=cpm(deg_norm,log = TRUE,prior.count = 2)

# visalization effect for normalization 

par(mfrow=c(1,2))

# befor normalization  

deg_cpm=cpm(deg_object_filter,log = TRUE,prior.count = 2) 

boxplot(deg_cpm,las=2,col=as.integer(plura_meta$condition),main="before TMMwsp normalization",ylab="log2_cpm",cex.axis=0.7) 

# after normalization 

boxplot(deg_norm_cpm,las=2,col=as.integer(plura_meta$condition),main="after TTMWSP normalization",ylab="log2_cpm",cex.axis=0.07)

# MDS PLOT

par(mfrow = c(1,1))

plotMDS(deg_norm_cpm,col=c("royalblue","firebrick")[as.integer(plura_meta$condition)],cex = 1.5,gene.selection = "common",main="MDS_PLOT")

legend("topright", legend = levels(plura_meta$condition),col    = c("royalblue", "firebrick"),pch    = 16)

#interactive MDS 

glimmaMDS(deg_norm_cpm,groups = plura_meta$condition)

# herarchical  cluster dendrogram

hc=hclust(dist(t(deg_norm_cpm)),method = "complete")

plot(hc,main="Hierarchical Clustering of Samples", sub  = "Euclidean distance, log2-CPM",xlab="",cex = 0.8)

#sample coleration heatmap

sample_colr=cor(deg_norm_cpm,method = "pearson")

pheatmap(sample_colr,color = colorRampPalette(c("white","blue"))(100),main = "sample_colleration heatmap",angle_col = 45,fontsize =10,las=2)

#denisty plot of log_cpm

df_long=melt(deg_norm_cpm,varnames = c("gene","sample"),value.name = "log_cpm")

df_long$condition=plura_meta[df_long$sample,"condition"]

ggplot(df_long,aes(x = log_cpm,colour = condition,group = sample))+geom_density(alpha = 0.5)+theme_bw()+labs(title = "log_cpm density per sample",x="logcpm")

# design matrix 

design=model.matrix(~0+condition,data = plura_meta)

colnames(design)

dim(design)

#check design rank 

qr(design)$rank==ncol(design)

pheatmap(design,cluster_rows = FALSE,cluster_cols = FALSE,main = "Design Matrix" ,color = c("white", "steelblue"))

#voom transformation 

voom_d=voom(deg_norm,design,plot = TRUE)

#voom with quailty weight

qualty_vw=voomWithQualityWeights(deg_norm,design = design,plot = TRUE)

qualty_vw=voomWithQualityWeights(deg_norm,design = design,plot = TRUE,col = as.integer(plura_meta$condition))

qualty_vw$weights

# linear model fitting

fit_lmm=lmFit(voom_d,design)

#apply embrical byas

fit_em=eBayes(fit_lmm)

cat("fit df prior:",fit_em$df.prior,"\n")

cat("prior variance:",fit_em$s2.prior,"\n")

# extract result 

conditionM_Epithelial=topTable(fit_em,coef = "conditionM_Epithelial",number = Inf,adjust.method = "BH") 

conditionM_mixed =topTable(fit_em,coef = "conditionM_mixed" ,number = Inf,adjust.method = "BH")

conditionM_Sarcomatoid =topTable(fit_em,coef = "conditionM_Sarcomatoid",number = Inf,adjust.method = "BH")

conditionMesothelioma=topTable(fit_em,coef = "conditionMesothelioma" ,number = Inf,adjust.method = "BH")

summary(decideTests(fit_em,p.value = 0.05))

#make constract

make_matrix=makeContrasts(epi_vs_mix= conditionM_Epithelial-conditionM_mixed,
                          epi_vs_Sar=conditionM_Epithelial-conditionM_Sarcomatoid,
                          epi_vs_Mes= conditionM_Epithelial-conditionMesothelioma,
                          mix_vs_Sar= conditionM_Sarcomatoid-conditionM_mixed ,
                          Sar_vs_Mes=conditionMesothelioma-conditionM_Sarcomatoid ,levels = design)

fit_matx=contrasts.fit(fit_lmm,make_matrix)

matx_em=eBayes(fit_matx)

#extract result

epi_vs_mix =topTable(matx_em,coef = "epi_vs_mix",number = Inf)

epi_vs_Sar=topTable(matx_em,coef = "epi_vs_Sar",number = Inf)

epi_vs_Mes=topTable(matx_em,coef = "epi_vs_Mes",number = Inf)

mix_vs_Sar=topTable(matx_em,coef = "mix_vs_Sar",number = Inf)

Sar_vs_Mes=topTable(matx_em,coef = "Sar_vs_Mes",number = Inf)

# add gene anotation epi_vs_mix

epi_vs_mix $ENTREZID =mapIds(org.Hs.eg.db,keys =rownames(epi_vs_mix ),column ="ENTREZID",keytype ="ENSEMBL",multiVals = "first")

epi_vs_mix $symbol =mapIds(org.Hs.eg.db,keys =rownames(epi_vs_mix ),column ="SYMBOL",keytype ="ENSEMBL",multiVals = "first")

epi_vs_mix $genename=mapIds(org.Hs.eg.db,keys =rownames(epi_vs_mix ),column ="GENENAME",keytype ="ENSEMBL",multiVals = "first")

# add gene anotation epi_vs_Sar 

epi_vs_Sar $ENTREZID =mapIds(org.Hs.eg.db,keys =rownames(epi_vs_Sar),column ="ENTREZID",keytype ="ENSEMBL",multiVals = "first")

epi_vs_Sar $symbol =mapIds(org.Hs.eg.db,keys =rownames(epi_vs_Sar),column ="SYMBOL",keytype ="ENSEMBL",multiVals = "first") 

epi_vs_Sar $genename =mapIds(org.Hs.eg.db,keys =rownames(epi_vs_Sar),column ="GENENAME",keytype ="ENSEMBL",multiVals = "first")

# add gene anotation epi_vs_Mes

epi_vs_Mes $ENTREZID =mapIds(org.Hs.eg.db,keys =rownames(epi_vs_Mes),column ="ENTREZID",keytype ="ENSEMBL",multiVals = "first")

epi_vs_Mes $symbol =mapIds(org.Hs.eg.db,keys =rownames(epi_vs_Mes),column ="SYMBOL",keytype ="ENSEMBL",multiVals = "first")

epi_vs_Mes $genename =mapIds(org.Hs.eg.db,keys =rownames(epi_vs_Mes),column ="GENENAME",keytype ="ENSEMBL",multiVals = "first")

# add gene anotation mix_vs_Sar 

mix_vs_Sar $ENTREZID =mapIds(org.Hs.eg.db,keys =rownames( mix_vs_Sar) ,column ="ENTREZID",keytype ="ENSEMBL",multiVals = "first")

mix_vs_Sar $symbol =mapIds(org.Hs.eg.db,keys =rownames( mix_vs_Sar),column ="SYMBOL",keytype ="ENSEMBL",multiVals = "first")

mix_vs_Sar $genename =mapIds(org.Hs.eg.db,keys =rownames( mix_vs_Sar),column ="GENENAME",keytype ="ENSEMBL",multiVals = "first")

# add gene anotation Sar_vs_Mes

Sar_vs_Mes$ENTREZID =mapIds(org.Hs.eg.db,keys =rownames(Sar_vs_Mes) ,column ="ENTREZID",keytype ="ENSEMBL",multiVals = "first")

Sar_vs_Mes$symbol =mapIds(org.Hs.eg.db,keys =rownames( Sar_vs_Mes),column ="SYMBOL",keytype ="ENSEMBL",multiVals = "first")

Sar_vs_Mes$genename=mapIds(org.Hs.eg.db,keys =rownames( Sar_vs_Mes),column ="GENENAME",keytype ="ENSEMBL",multiVals = "first")

# delete NA 
add_annotation <- function(df) {
  df$ENTREZID <- mapIds(org.Hs.eg.db,
                        keys = rownames(df),
                        column = "ENTREZID",
                        keytype = "ENSEMBL",
                        multiVals = "first")
  
  df$symbol <- mapIds(org.Hs.eg.db,
                      keys = rownames(df),
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")
  
  df$genename <- mapIds(org.Hs.eg.db,
                        keys = rownames(df),
                        column = "GENENAME",
                        keytype = "ENSEMBL",
                        multiVals = "first")
  
  df <- df[complete.cases(df[, c("ENTREZID", "symbol", "genename")]), ]
  
  return(df)
}

epi_vs_mix <- add_annotation(epi_vs_mix)
epi_vs_Sar <- add_annotation(epi_vs_Sar)
epi_vs_Mes <- add_annotation(epi_vs_Mes)
mix_vs_Sar <- add_annotation(mix_vs_Sar)
Sar_vs_Mes <- add_annotation(Sar_vs_Mes)

#save results

write.csv(epi_vs_mix,"epi_vs_mix.csv",row.names = TRUE)

write.csv(epi_vs_Sar,"epi_vs_Sar.csv",row.names = TRUE)

write.csv(epi_vs_Mes,"epi_vs_Mes.csv",row.names = TRUE)

write.csv(Sar_vs_Mes,"Sar_vs_Mes.csv",row.names = TRUE)

write.csv(mix_vs_Sar,"mix_vs_Sar.csv",row.names = TRUE)

# significant gene

sig_epi_vs_Sar=epi_vs_Sar  %>% filter(adj.P.Val < 0.05, abs(logFC) > 1)

sig_epi_vs_Mes=epi_vs_Mes %>% filter(adj.P.Val < 0.05, abs(logFC) >1)

sig_mix_vs_Sar=mix_vs_Sar %>% filter(adj.P.Val < 0.05, abs(logFC) >1)

# significant  up & down gene sig_sar_ep

up_sig_mix_vs_Sar=filter(sig_mix_vs_Sar,logFC >1)

down_sig_mix_vs_Sar=filter(sig_mix_vs_Sar,logFC < 1)

up_sig_epi_vs_Sar=filter(sig_epi_vs_Sar,logFC > 1)

up_sig_epi_vs_Mes=filter(sig_epi_vs_Mes,logFC > 1)

down_sig_epi_vs_Mes=filter(sig_epi_vs_Mes,logFC < 1)


#  save significant gene


write.csv(sig_epi_vs_Sar,"sig_epi_vs_Sar.csv",row.names = TRUE)

write.csv(sig_epi_vs_Mes,"sig_epi_vs_Mes.csv",row.names = TRUE)

write.csv(sig_mix_vs_Sar,"sig_mix_vs_Sar.csv",row.names = TRUE)

#  save significant  up & down gene sig_sar_ep

write.csv(up_sig_mix_vs_Sar,"up_sig_mix_vs_Sar.csv",row.names = TRUE)

write.csv(down_sig_mix_vs_Sar,"down_sig_mix_vs_Sar.csv",row.names = TRUE)

write.csv(up_sig_epi_vs_Sar,"up_sig_epi_vs_Sar.csv",row.names = TRUE)

write.csv(up_sig_epi_vs_Mes,"up_sig_epi_vs_Mes.csv",row.names = TRUE)

write.csv(down_sig_epi_vs_Mes,"down_sig_epi_vs_Mes.csv",row.names = TRUE)

#VISUALIZATION
v2 <- voom(deg_norm, design, plot = TRUE,  span = 0.5)

# MA PLOT

plotMD(fit_matx,
              column = "epi_vs_mix",
             status = epi_vs_mix$adj.P.Val < 0.05,
             values = c(1, 0),
             col    = c("red", "black"),
              main   = "MA Plot - epi_vs_mix",
             ylab   = "log2 Fold Change") 
            abline(h = 0, col = "blue", lty = 2)
           
            
plotMD(fit_matx,
                  column = "epi_vs_Sar",
                   status =epi_vs_Sar$adj.P.Val < 0.05,
                   values = c(1, 0),
                   col    = c("red", "black"),
                   main   = "MA Plot - epi_vs_Sar",
                   ylab   = "log2 Fold Change") 
            abline(h = 0, col = "blue", lty = 2)        


plotMD(fit_matx,
                   column = "epi_vs_Mes",
                   status =epi_vs_Mes$adj.P.Val < 0.05,
                   values = c(1, 0),
                   col    = c("red", "black"),
                   main   = "MA Plot - epi_vs_Sar",
                   ylab   = "log2 Fold Change") 
            abline(h = 0, col = "blue", lty = 2)        
            


            
plotMD(fit_matx,
                  column = "mix_vs_Sar",
                   status =mix_vs_Sar$adj.P.Val < 0.05,
                   values = c(1, 0),
                   col    = c("red", "black"),
                   main   = "MA Plot - mix_vs_Sar",
                   ylab   = "log2 Fold Change") 
            abline(h = 0, col = "blue", lty = 2)    


            
plotMD(fit_matx,
                column = "Sar_vs_Mes",
                status =mix_vs_Sar$adj.P.Val < 0.05,
                values = c(1, 0),
                col    = c("red", "black"),
                main   = "MA Plot - Sar_vs_Mes",
                ylab   = "log2 Fold Change") 
            abline(h = 0, col = "blue", lty = 2)    



# volcano plot 

 EnhancedVolcano(epi_vs_mix,
                            lab      = epi_vs_mix$symbol,
                            x        = "logFC",
                            y        = "adj.P.Val",
                            title    = "epi_vs_mix",
                            subtitle = "Limma-voom | epi_vs_mix ",
                            pCutoff  = 0.05,
                            FCcutoff = 1,
                            col      = c("grey40", "forestgreen", "royalblue", "red2")
            )

 
 EnhancedVolcano(epi_vs_Sar,
                 lab      = epi_vs_Sar$symbol,
                 x        = "logFC",
                 y        = "adj.P.Val",
                 title    = "epi_vs_Sar",
                 subtitle = "Limma-voom | epi_vs_Sar",
                 pCutoff  = 0.05,
                 FCcutoff = 1,
                 col      = c("grey40", "forestgreen", "royalblue", "red2")
 )

 
 EnhancedVolcano(epi_vs_Mes,
                 lab      = epi_vs_Mes$symbol,
                 x        = "logFC",
                 y        = "adj.P.Val",
                 title    = "epi_vs_Mes",
                 subtitle = "Limma-voom | epi_vs_Mes",
                 pCutoff  = 0.05,
                 FCcutoff = 1,
                 col      = c("grey40", "forestgreen", "royalblue", "red2")
 )
 
 
 EnhancedVolcano(mix_vs_Sar,
                 lab      = mix_vs_Sar$symbol,
                 x        = "logFC",
                 y        = "adj.P.Val",
                 title    = "epi_vs_Mes",
                 subtitle = "Limma-voom | mix_vs_Sar",
                 pCutoff  = 0.05,
                 FCcutoff = 1,
                 col      = c("grey40", "forestgreen", "royalblue", "red2")
 )
 
 EnhancedVolcano(Sar_vs_Mes,
                 lab      =Sar_vs_Mes$symbol,
                 x        = "logFC",
                 y        = "adj.P.Val",
                 title    = "epi_vs_Mes",
                 subtitle = "Limma-voom | Sar_vs_Mes",
                 pCutoff  = 0.05,
                 FCcutoff = 1,
                 col      = c("grey40", "forestgreen", "royalblue", "red2")
 )
 
 
# heatmap sig_epi_vs_Sar 
 
top_epi_vs_Sar=rownames(head(sig_epi_vs_Sar,50))

mat_epi_vs_Sar=v2$E[top_epi_vs_Sar,]

mat_epi_vs_Sar=mat_epi_vs_Sar-rowMeans(mat_epi_vs_Sar)
 
rownames(mat_epi_vs_Sar)= sig_epi_vs_Sar[top_epi_vs_Sar,"symbol"]
 
pheatmap(mat_epi_vs_Sar,
         annotation_col = data.frame(
          Condition = plura_meta$condition,
           row.names = colnames(mat_epi_vs_Sar)
         ),
         color    = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         scale    = "row",
         main     = "Top 50 DE Genes (Limma-voom epi_vs_Sar )",
         fontsize_row = 7
)


# heatmap sig_epi_vs_Mes


top_sig_epi_vs_Mes=rownames(head(sig_epi_vs_Mes,50))

mat_sig_epi_vs_Mes=v2$E[top_sig_epi_vs_Mes,]

mat_sig_epi_vs_Mes=mat_sig_epi_vs_Mes-rowMeans(mat_sig_epi_vs_Mes)

rownames(mat_sig_epi_vs_Mes)= sig_epi_vs_Mes[top_sig_epi_vs_Mes,"symbol"]

pheatmap(mat_sig_epi_vs_Mes,
         annotation_col = data.frame(
           Condition = plura_meta$condition,
           row.names = colnames(mat_sig_epi_vs_Mes)
         ),
         color    = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         scale    = "row",
         main     = "Top 50 DE Genes (Limma-voom epi_vs_Mes )",
         fontsize_row = 7
)


# heatmap sig_mix_vs_Sar



top_sig_mix_vs_Sar=rownames(head(sig_mix_vs_Sar,50))

mat_sig_mix_vs_Sar=v2$E[top_sig_mix_vs_Sar,]

mat_sig_mix_vs_Sar=mat_sig_mix_vs_Sar-rowMeans(mat_sig_mix_vs_Sar)

rownames(mat_sig_mix_vs_Sar)=sig_mix_vs_Sar[top_sig_mix_vs_Sar,"symbol"]

pheatmap(mat_sig_mix_vs_Sar,
         annotation_col = data.frame(
           Condition = plura_meta$condition,
           row.names = colnames(mat_sig_mix_vs_Sar)
         ),
         color    = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         scale    = "row",
         main     = "Top 50 DE Genes (Limma-voom sig_mix_vs_Sar )",
         fontsize_row = 7
)

#venn Diagram of overlapping DE genes 

result_mat=decideTests(matx_em,
                        method = "separate",
                       adjust.method = "BH",
                       p.value  = 0.05,
                       lfc      = 1  )
vennDiagram(result_mat,
            include = c("up", "down"),
            counts.col = c("red", "blue"),
            circle.col = c("#1B5E20", "#1A237E", "#6A1B9A"),
            cex = c(1, 1.2, 1)
)


# GENE set enrichment analysis epi_vs_Sar

gene_epi_vs_Sar_list=epi_vs_Sar$t

names(gene_epi_vs_Sar_list)=epi_vs_Sar$ENTREZID 

gene_epi_vs_Sar_list = gene_epi_vs_Sar_list[!is.na(names(gene_epi_vs_Sar_list))]

gene_epi_vs_Sar_list= sort(gene_epi_vs_Sar_list, decreasing = TRUE)

 gsea_epi_vs_Sar <- gseKEGG(
  geneList = gene_epi_vs_Sar_list,
  organism = "hsa",
  minGSSize = 15,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

 #  GENE set enrichment analysis epi_vs_mix
 
 gene_epi_vs_mix_list=epi_vs_mix$t
 
 names(gene_epi_vs_mix_list)=epi_vs_mix$ENTREZID 
 
 gene_epi_vs_mix_list= gene_epi_vs_mix_list[!is.na(names(gene_epi_vs_mix_list))]
 
 gene_epi_vs_mix_list= sort( gene_epi_vs_mix_list, decreasing = TRUE)
 
 gsea_epi_vs_mix <- gseKEGG(
   geneList =gene_epi_vs_mix_list,
   organism = "hsa",
   minGSSize = 15,
   maxGSSize = 500,
   pvalueCutoff = 0.05,
   verbose = FALSE
 )
 
 # #  GENE set enrichment analysis epi_vs_Mes
 
 gene_epi_vs_Mes_list=epi_vs_Mes$t
 
 names( gene_epi_vs_Mes_list)=epi_vs_Mes$ENTREZID 
 
 gene_epi_vs_Mes_list = gene_epi_vs_Mes_list[!is.na(names(gene_epi_vs_Mes_list))]
 
 gene_epi_vs_Mes_list= sort( gene_epi_vs_Mes_list, decreasing = TRUE)
 
 gsea_epi_vs_Mes <- gseKEGG(
   geneList =  gene_epi_vs_Mes_list,
   organism = "hsa",
   minGSSize = 15,
   maxGSSize = 500,
   pvalueCutoff = 0.05,
   verbose = FALSE
 )
 
 # #  GENE set enrichment analysis mix_vs_Sar
 
 gene_mix_vs_Sar_list= mix_vs_Sar$t
 
 names(  gene_mix_vs_Sar_list)=mix_vs_Sar$ENTREZID 
 
 gene_mix_vs_Sar_list = gene_mix_vs_Sar_list[!is.na(names(gene_mix_vs_Sar_list))]
 
 gene_mix_vs_Sar_list= sort(  gene_mix_vs_Sar_list, decreasing = TRUE)
 
 gsea__mix_vs_Sar <- gseKEGG(
   geneList =  gene_mix_vs_Sar_list,
   organism = "hsa",
   minGSSize = 15,
   maxGSSize = 500,
   pvalueCutoff = 0.05,
   verbose = FALSE
 )
 
 #  GENE set enrichment analysis  Sar_vs_Mes
 gene_Sar_vs_Mes_list=  Sar_vs_Mes$t
 
 names(   gene_Sar_vs_Mes_list)=Sar_vs_Mes$ENTREZID 
 
 gene_Sar_vs_Mes_list = gene_Sar_vs_Mes_list[!is.na(names(gene_Sar_vs_Mes_list))]
 
 gene_Sar_vs_Mes_list= sort( gene_Sar_vs_Mes_list, decreasing = TRUE)
 
 gsea__Sar_vs_Mes <- gseKEGG(
   geneList =  gene_Sar_vs_Mes_list,
   organism = "hsa",
   minGSSize = 15,
   maxGSSize = 500,
   pvalueCutoff = 0.05,
   verbose = FALSE
 )
 
 
 # save gene set enrichment analysis
 

  write.csv( gsea_epi_vs_Sar," gsea_epi_vs_Sar.csv",row.names = FALSE)

 write.csv(  gsea_epi_vs_mix,"  gsea_epi_vs_mix.csv",row.names = FALSE)
 
 write.csv(   gsea_epi_vs_Mes ,"   gsea_epi_vs_Mes .csv",row.names = FALSE)
 
 write.csv(    gsea__Sar_vs_Mes ,"gsea__Sar_vs_Mes.csv",row.names = FALSE)
 
 write.csv(    gsea__mix_vs_Sar ," gsea__mix_vs_Sar.csv",row.names = FALSE)
 
 # visualiztion Sar_vs_Mes
 
 dotplot(gsea__Sar_vs_Mes , showCategory = 15)
 
 ridgeplot(gsea__Sar_vs_Mes ,showCategory = 15)
 
 emapplot(pairwise_termsim(gsea__Sar_vs_Mes) ,,showCategory = 15)
 
 
 # visualiztion gsea__mix_vs_Sar
 
 
 dotplot(gsea__mix_vs_Sar , showCategory = 15)
 
 ridgeplot(gsea__mix_vs_Sar,showCategory = 15)
 
 emapplot(pairwise_termsim(gsea__mix_vs_Sar) ,,showCategory = 20)
 
 # visualiztion  gsea_epi_vs_Mes 
 
 dotplot(gsea_epi_vs_Mes  , showCategory = 15)
 
 ridgeplot(gsea_epi_vs_Mes ,showCategory = 15)
 
 emapplot(pairwise_termsim(gsea_epi_vs_Mes ) ,,showCategory = 20)
 
# visualiztion  gsea_epi_vs_mix  
 
 dotplot( gsea_epi_vs_mix  , showCategory = 15)
 
 ridgeplot( gsea_epi_vs_mix,showCategory = 15)
 
 emapplot(pairwise_termsim( gsea_epi_vs_mix ) ,,showCategory = 20)
 
 # visualiztion  epi_vs_Sar
 
 
 dotplot(  gsea_epi_vs_Sar , showCategory = 15)
 
 ridgeplot(  gsea_epi_vs_Sar,showCategory = 15)
 
 emapplot(pairwise_termsim(  gsea_epi_vs_Sar ) ,,showCategory = 20)
 
 # reactome sig_epi_vs_Sar
 
 reactome_epi_vs_Sar=enrichPathway(gene =unique(na.omit(sig_epi_vs_Sar$ENTREZID)),organism = "human",pvalueCutoff = 0.05)
 
 # reactome sig_epi_vs_Sar
 
 reactome_epi_vs_Mes=enrichPathway(gene =unique(na.omit(sig_epi_vs_Mes$ENTREZID)),organism = "human",pvalueCutoff = 0.05)
 
 # reactome sig_epi_vs_Sar
 
 reactome_sig_mix_vs_Sar =enrichPathway(gene =unique(na.omit(sig_mix_vs_Sar $ENTREZID)),organism = "human",pvalueCutoff = 0.05)
 
 
 # save reactome result
 
 write.csv( reactome_epi_vs_Sar," reactome_epi_vs_Sar.csv",row.names = FALSE)
 



 
 