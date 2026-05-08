# Load packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(IHW)
library(GSVA)
library(car)
library(dunn.test)

# Aggregate counts
Idents(macs)<-'disease'
macs<-subset(macs, idents=c('Emphysema', 'Fibrosis', 'Tu.free'))

table(macs$integrated_snn_res.0.2, macs$LuT)

Idents(macs)<-macs$integrated_snn_res.0.2
final_list<-list()

for(i in levels(macs$integrated_snn_res.0.2)){
  
  tmp<-subset(macs, idents=i)
  keep_lut<-names(which(table(tmp$LuT) >= 3))
  if(length(keep_lut) < 1) next
  
  tmp<-subset(tmp, subset = LuT %in% keep_lut)
  
  pseudobulk<-AggregateExpression(tmp,
                                  group.by="LuT",
                                  assays="RNA",
                                  slot="counts")
  
  pseudobulk<-pseudobulk$RNA
  pseudobulk<-as.data.frame(pseudobulk)
  
  metadata<-data.frame(patient=colnames(pseudobulk))
  metadata$disease<-'Fibrosis'
  idx<-metadata$patient %in% c('LuT_314','LuT_349','LuT_383','LuT_410','LuT_433','LuT_736')
  if(any(idx)) metadata$disease[idx]<-'Emphysema'
  idx<-metadata$patient %in% c('LuT_437','LuT_553','LuT_580')
  if(any(idx)) metadata$disease[idx]<-'Control'
  rownames(metadata)<-colnames(pseudobulk)
  
  if(any(table(factor(metadata$disease, levels=c("Emphysema", "Fibrosis", "Control"))) < 3)) next
  
  dds<-DESeqDataSetFromMatrix(countData=pseudobulk,
                              colData=metadata,
                              design=~disease)
  
  keep<-rowSums(counts(dds) >= 10) >= 3
  dds<-dds[keep,]
  dds<-DESeq(dds)
  
  res.emph<-results(dds,
                    contrast=list("disease_Emphysema_vs_Control","disease_Fibrosis_vs_Control"),
                    listValues=c(1, -0.5),
                    lfcThreshold=0,
                    alpha=0.1,
                    filterFun=ihw,
                    altHypothesis="greaterAbs")
  
  res.emph<-as.data.frame(res.emph)
  res.emph<-res.emph[res.emph$pvalue <= 0.05, ]
  res.emph<-res.emph[!is.na(res.emph$baseMean),]
  
  if(dim(res.emph)[1]>0){
    res.emph$gene<-rownames(res.emph)
    res.emph$direction<-'Emphysema vs rest'
    res.emph$cluster<-paste0('macs', i)
  }
  
  res.fib<-results(dds,
                   contrast=list("disease_Fibrosis_vs_Control","disease_Emphysema_vs_Control"),
                   listValues=c(1, -0.5),
                   lfcThreshold=0,
                   alpha=0.1,
                   filterFun=ihw,
                   altHypothesis="greaterAbs")
  
  res.fib<-as.data.frame(res.fib)
  res.fib<-res.fib[res.fib$pvalue <= 0.05, ]
  res.fib<-res.fib[!is.na(res.fib$baseMean),]
  
  if(dim(res.fib)[1]>0){
    res.fib$gene<-rownames(res.fib)
    res.fib$direction<-'Fibrosis vs rest'
    res.fib$cluster<-paste0('macs', i)
  }
  
  tmp2<-rbind(res.emph, res.fib)
  final_list[[i]]<-tmp2
}

df<-do.call(rbind, final_list)

write.csv(df, 'macs_pare.csv')