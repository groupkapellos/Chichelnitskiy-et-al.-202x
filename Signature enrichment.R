# Load packages
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(GSVA)
library(reshape2)
library(ggplot2)
library(car)
library(dunn.test)
library(ggvenn)

# Load data
gse37768<-read.delim('GSE37768.txt')
rownames(gse37768)<-gse37768$Symbol
gse37768$Symbol<-NULL

metadata<-read.csv('GSE37768_metadata.csv')
rownames(metadata)<-metadata$geo_accession

gse48149<-read.delim('GSE48149.txt')
rownames(gse48149)<-gse48149$Symbol
gse48149$Symbol<-NULL

metadata<-read.csv('GSE48149_metadata.csv')
rownames(metadata)<-metadata$geo_accession
metadata$X.Sample_source_name_ch1<-gsub('SSc PF', 'SSc-PF', metadata$X.Sample_source_name_ch1)

# Select gene signatures
setwd('D:/Theo/Documents/LIMES/LIMES 2019-20+/Papers/Hannover paper/Tables')
#genes<-read_xlsx('Table S3.xlsx', sheet=4)
#genes<-read_xlsx('Table S3.xlsx', sheet=5)
#genes<-read_xlsx('Table S3.xlsx', sheet=6)
#genes<-read_xlsx('Table S4.xlsx', sheet=3)
#genes<-read_xlsx('Table S5.xlsx', sheet=4)
#genes<-read_xlsx('Table S5.xlsx', sheet=5)
genes<-read_xlsx('Table S5.xlsx', sheet=6)

genelist<-list()
for (i in unique(genes$Cluster)){
  genelist[[i]]<-unique(genes[genes$Cluster==i & genes$Disease=='Emphysema',]$Gene)
}

#genelist<-list()
#for (i in unique(genes$Cluster)){
#  genelist[[i]]<-unique(genes[genes$Cluster==i & genes$Disease=='Fibrosis',]$Gene)
#}

# Run GSVA
es<-gsva(as.matrix(gse37768),
         genelist,
         min.sz=5,
         max.sz=500,
         verbose=T,
         method='ssgsea',
         mx.diff=F,               
         parallel.sz=1)

es<-gsva(as.matrix(gse48149),
         genelist,
         min.sz=5,
         max.sz=500,
         verbose=T,
         method='ssgsea',
         mx.diff=F,               
         parallel.sz=1)

es<-melt(es)
es$stage<-metadata$Stage[match(es$Var2, metadata$geo_accession)]
es$stage<-factor(es$stage, levels=c('nonsmoker', 'smoker', 'COPD'))

es<-melt(es)
es$stage<-metadata$X.Sample_source_name_ch1[match(es$Var2, metadata$X.Sample_geo_accession)]
es$stage<-factor(es$stage, levels=c('NL','IPF','SSc-PF','IPAH (PPH)','SSc-PAH'))

ggplot(es, aes(x=stage, y=value)) +
  geom_boxplot() +
  facet_wrap(~Var1, scale='free') +
  theme_bw()

for(i in unique(es$Var1)){
  print(i)
  tmp<-es[es$Var1==i,]
  res.aov<-aov(value~stage, tmp)
  var<-leveneTest(value~stage, tmp)
  var<-var$`Pr(>F)`[1]
  resi<-shapiro.test(residuals(object=res.aov))
  resi<-resi$p.value
  
  if(var | resi<=0.05){
    print('Kruskal-Wallis')
    print(oneway.test(value ~ stage, data=tmp))
    print(pairwise.t.test(tmp$value, tmp$stage, p.adjust.method="BH", pool.sd=FALSE))
    
  }else{
    print('ANOVA')
    print(summary(res.aov))
    print(TukeyHSD(res.aov))
  }
}
