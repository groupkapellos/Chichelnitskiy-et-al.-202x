# Load packages
library(Seurat)
library(dplyr)
library(org.Hs.eg.db)
library(biomaRt)
library(clusterProfiler)
library(ggplot2)

## Parenchyma
# Compare immune cell populations between emphysema-fibrosis
# alveolar macs: 0,1,2,4
# alveolar macs: 5 (no)
# interstitial macs: 3
# monos: 0,1,2,3,4
# nonclassical monos: 5
# cdcs: 0,1
# pdcs: 2
# dcs: 3,4 (no)
# neutros: 0,1
# neutros: 4,5 (no)
# imm neutros: 2,3
# mast: 1,2
# mast: 0 (no)
# cyto: 0
# trm: 1
# trm-like: 2,3,6
# tcm: 5
# cycling: 7
# trm-like: 4 (no)
# tem: 8,9 (no)
# activated: 10 (no)
# b mature: 1
# plasma igg: 2
# plasma iga: 3
# b: 0,4 (no)
# nks: 0,1,2,3

Idents(monos)<-'integrated_snn_res.0.3'

deg2<-data.frame()
for(i in c(0)){
  test<-subset(monos, idents=i)
  Idents(test)<-'disease'
  deg<-FindMarkers(test, ident.1='Emphysema', ident.2='Fibrosis', min.pct=0.2, logfc.threshold=0.25, test.use='MAST', only.pos=T)
  deg<-deg %>% arrange(-avg_log2FC)
  deg$cluster<-paste0('Cluster_',i)
  deg$disease<-'Emphysema'
  deg$gene<-rownames(deg)
  deg2<-rbind(deg2,deg)
  
  deg<-FindMarkers(test, ident.1='Fibrosis', ident.2='Emphysema', min.pct=0.2, logfc.threshold=0.25, test.use='MAST', only.pos=T)
  deg<-deg %>% arrange(-avg_log2FC)
  deg$cluster<-paste0('Cluster_',i)
  deg$disease<-'Fibrosis'
  deg$gene<-rownames(deg)
  deg2<-rbind(deg2,deg)
}
deg2<-deg2[deg2$p_val_adj<=0.05,]
write.csv(deg2, 'monos.csv')

GO.up<-enrichGO(gene=data[data$Cluster=='2' & data$Disease=='Fibrosis',]$Gene,
                OrgDb=org.Hs.eg.db,
                ont='BP',
                keyType='SYMBOL',
                minGSSize=10,
                maxGSSize=500,
                pvalueCutoff=0.05,
                readable=FALSE)

GO.up<-simplify(GO.up, cutoff=0.7, by="p.adjust", select_fun=min)
head(GO.up)

df<-as.data.frame(GO.up)
df<-df[df$Count>=3,]
df$cluster<-'Cluster 2'
df$disease<-'Fibrosis'

#deg2<-data.frame()
deg2<-rbind(deg2, df)

## Lymph nodes
# Compare disease status among immune cell populations
# macs: 0,1
# monos: 0
# cdcs: 0
# pdcs: 1
# resting: 0
# activated: 1,2
# temra: 3
# memory: 4,5
# cycling: 6
# b: 1
# plasma: 0,2
# nk: 1

Idents(macs)<-'integrated_snn_res.0.3'

deg2<-data.frame()
for(i in c(0,1)){
  test<-subset(macs, idents=i)
  Idents(test)<-'disease'
  deg<-FindMarkers(test, ident.1='Emphysema', ident.2='Fibrosis', min.pct=0.2, logfc.threshold=0.25, test.use='MAST', only.pos=T)
  deg<-deg %>% arrange(-avg_log2FC)
  deg$cluster<-paste0('Cluster_',i)
  deg$disease<-'Emphysema'
  deg$gene<-rownames(deg)
  deg2<-rbind(deg2,deg)
  
  deg<-FindMarkers(test, ident.1='Fibrosis', ident.2='Emphysema', min.pct=0.2, logfc.threshold=0.25, test.use='MAST', only.pos=T)
  deg<-deg %>% arrange(-avg_log2FC)
  deg$cluster<-paste0('Cluster_',i)
  deg$disease<-'Fibrosis'
  deg$gene<-rownames(deg)
  deg2<-rbind(deg2,deg)
}
deg2<-deg2[deg2$p_val_adj<=0.05,]
write.csv(deg2, 'nks.csv')

data<-deg2
GO.up<-enrichGO(gene=data[data$cluster=='Cluster_0' & data$disease=='Emphysema',]$gene,
                OrgDb=org.Hs.eg.db,
                ont='BP',
                keyType='SYMBOL',
                minGSSize=10,
                maxGSSize=500,
                pvalueCutoff=0.05,
                readable=FALSE)

GO.up<-simplify(GO.up, cutoff=0.7, by="p.adjust", select_fun=min)
head(GO.up)

df<-as.data.frame(GO.up)
df<-df[df$Count>=3,]
df$cluster<-'Cluster 1'
df$disease<-'Fibrosis'

#deg2<-data.frame()
deg2<-rbind(deg2, df)
