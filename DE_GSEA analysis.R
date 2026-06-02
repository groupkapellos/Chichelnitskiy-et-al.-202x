# Load packages
library(Seurat)
library(readxl)
library(dplyr)
library(org.Hs.eg.db)
library(biomaRt)
library(clusterProfiler)
library(ggplot2)

# Set directory
setwd('C:/Users/theodoros.kapellos/Downloads')
getwd()

# Load DE genes
genes<-readxl::read_xlsx('Table S5.xlsx', sheet=6)

# Perform pathway analysis
deg2<-data.frame()

for(i in unique(genes$cluster)[c(1)]){
  tmp1<-genes[genes$direction=='Higher in emphysema' & genes$cluster==i,]$Gene
  
  GO.up<-enrichGO(gene=tmp1,
                  OrgDb=org.Hs.eg.db,
                  ont='BP',
                  keyType='SYMBOL',
                  minGSSize=10,
                  maxGSSize=500,
                  pvalueCutoff=0.05,
                  readable=FALSE)
  
  GO.up<-simplify(GO.up, cutoff=0.7, by="p.adjust", select_fun=min)
  df1<-as.data.frame(GO.up)
  
  if (nrow(df1) > 0) {
    df1<-df1[df1$Count>=3,]
    df1$cluster<-i
    df1$disease<-'Higher in emphysema'
    }
  
  tmp2<-genes[genes$direction=='Higher in fibrosis' & genes$cluster==i,]$Gene
  
  GO.up2<-enrichGO(gene=tmp2,
                  OrgDb=org.Hs.eg.db,
                  ont='BP',
                  keyType='SYMBOL',
                  minGSSize=10,
                  maxGSSize=500,
                  pvalueCutoff=0.05,
                  readable=FALSE)
  
  GO.up2<-simplify(GO.up, cutoff=0.7, by="p.adjust", select_fun=min)
  df2<-as.data.frame(GO.up2)
  
  if (nrow(df2) > 0) {
    df2<-df2[df2$Count>=3,]
    df2$cluster<-i
    df2$disease<-'Higher in fibrosis'
  }
  
  deg<-rbind(df1, df2)
  deg2<-rbind(deg2, deg)
}

deg2<-deg2 %>% arrange(cluster, disease, -Count)
deg_final<-deg2 %>% group_by(cluster, disease) %>% arrange(desc(Count)) %>% slice_head(n=5)
