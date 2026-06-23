# Load packages
library(readxl)
library(GSVA)
library(reshape2)
library(ggplot2)
library(car)
library(dunn.test)
library(dplyr)
library(tidytext)

# Select gene signatures
#genes<-read_xlsx('Table S3.xlsx', sheet=4)
#genes<-read_xlsx('Table S3.xlsx', sheet=5)
#genes<-read_xlsx('Table S3.xlsx', sheet=6)
#genes<-read_xlsx('Table S4.xlsx', sheet=3)
#genes<-read_xlsx('Table S5.xlsx', sheet=4)
#genes<-read_xlsx('Table S5.xlsx', sheet=5)
genes<-read_xlsx('Table S5.xlsx', sheet=6)

genelist<-list()
for (i in unique(genes$cluster)){
  genelist[[i]]<-unique(genes[genes$cluster==i & genes$direction=='Higher in fibrosis',]$Gene)
}

# Load COPD data
gse37768<-read.delim('GSE37768.txt')
rownames(gse37768)<-gse37768$Symbol
gse37768$Symbol<-NULL

metadata<-read.csv('GSE37768_metadata.csv')
rownames(metadata)<-metadata$geo_accession

# Run GSVA
es<-gsva(as.matrix(gse37768),
         genelist,
         min.sz=5,
         max.sz=600,
         verbose=T,
         method='ssgsea',
         mx.diff=F,               
         parallel.sz=1)

es<-melt(es)
es$stage<-metadata$Stage[match(es$Var2, metadata$geo_accession)]
es$stage<-factor(es$stage, levels=c('nonsmoker', 'smoker', 'COPD'))

# Load IPF data
gse48149<-read.delim('GSE48149.txt')
rownames(gse48149)<-gse48149$Symbol
gse48149$Symbol<-NULL

metadata<-read.csv('GSE48149_metadata.csv')
rownames(metadata)<-metadata$geo_accession
metadata$X.Sample_source_name_ch1<-gsub('SSc PF', 'SSc-PF', metadata$X.Sample_source_name_ch1)

# Run GSVA
es2<-gsva(as.matrix(gse48149),
          genelist,
          min.sz=5,
          max.sz=600,
          verbose=T,
          method='ssgsea',
          mx.diff=F,       
          parallel.sz=1)

es2<-melt(es2)
es2$stage<-metadata$X.Sample_source_name_ch1[match(es2$Var2, metadata$X.Sample_geo_accession)]
es2$stage<-factor(es2$stage, levels=c('NL','IPF','SSc-PF','IPAH (PPH)','SSc-PAH'))

final<-rbind (es, es2)
final<-final[final$stage %in% c('COPD','IPF'),]

#results<-data.frame()

for(i in unique(final$Var1)){
  tmp<-final[final$Var1==i,]
  
  groups<-unique(tmp$stage)
  norm_p<-rep(NA, 2)
  names(norm_p)<-groups
  
  for (g in unique(tmp$stage)) {
    vals<-tmp[tmp$stage == g, "value"]
    vals<-as.data.frame(vals)
    norm_p[g]<-shapiro.test(vals$vals)$p.value
  }
  normal<-all(norm_p > 0.05)

  var<-leveneTest(value ~ stage, data=final)$`Pr(>F)`[1]
  
  if (normal) {
    
    if (var>0.05) {
      stat_p<-t.test(value ~ stage, data=tmp, var.equal=TRUE, alternative='two.sided')
      tmp2<-data.frame(parameter=i, meanCOPD=stat_p$estimate[1], meanIPF=stat_p$estimate[2], de='Fibrosis DE genes', test=stat_p$method, p.value=stat_p$p.value)
      results<-rbind(results, tmp2)
      
    } else {
      stat_p<-t.test(value ~ stage, data=tmp, var.equal=FALSE, alternative='two.sided')
      tmp2<-data.frame(parameter=i, meanCOPD=stat_p$estimate[1], meanIPF=stat_p$estimate[2], de='Fibrosis DE genes', test=stat_p$method, p.value=stat_p$p.value)
      results<-rbind(results, tmp2)
    }
    
  } else {
    stat_p<-wilcox.test(value ~ stage, data=tmp, alternative='two.sided')
    tmp2<-data.frame(parameter=i, meanCOPD=mean(tmp[tmp$stage == "COPD", "value"]), meanIPF=mean(tmp[tmp$stage == "IPF", "value"]), de='Fibrosis DE genes', test=stat_p$method, p.value=stat_p$p.value)
    results<-rbind(results, tmp2)
  }
}

results$direction<-ifelse(results$meanCOPD > results$meanIPF, 'Higher in COPD', 'Higher in IPF')
results$foldchange<-log2(results$meanIPF / results$meanCOPD)
results<-results %>% group_by(de) %>% mutate(adj.p.value=p.adjust(p.value, method="BH")) %>% ungroup()
results<-results[results$adj.p.value <= 0.05,]
