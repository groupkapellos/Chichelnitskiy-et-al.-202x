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

# Load datasets
Idents(combined_list)<-'LuT'
table(combined_list$LuT)

combined_list<-subset(combined_list, idents=unique(combined_list$LuT)[c(1:3, 5:16)])

# Aggregate counts
pseudobulk<-AggregateExpression(combined_list,
                                group.by="LuT",
                                assays="RNA",
                                slot="counts")

pseudobulk<-pseudobulk$RNA
pseudobulk<-as.data.frame(pseudobulk)

metadata<-data.frame(colnames(pseudobulk))

metadata$group<-'COPD'
metadata[metadata$colnames.pseudobulk. %in% c('LuT_437','LuT_553','LuT_580'),]$group<-'Control'
metadata[metadata$colnames.pseudobulk. %in% c('LuT_264','LuT_247','LuT_421','LuT_500','LuT_266','LuT_294','LuT_264'),]$group<-'IPF'
metadata[metadata$colnames.pseudobulk.%in% c('LuT_314','LuT_736'),]$group<-'CPFE'

metadata$disease<-'COPD'
metadata[metadata$colnames.pseudobulk. %in% c('LuT_437','LuT_553','LuT_580'),]$disease<-'Control'
metadata[metadata$colnames.pseudobulk. =='LuT_264',]$disease<-'Sjogren'
metadata[metadata$colnames.pseudobulk. %in% c('LuT_349','LuT_383'),]$disease<-'A1AT'
metadata[metadata$colnames.pseudobulk. %in% c('LuT_247','LuT_421','LuT_500'),]$disease<-'IPF'
metadata[metadata$colnames.pseudobulk. %in% c('LuT_266','LuT_294'),]$disease<-'EAA'
metadata[metadata$colnames.pseudobulk.%in% c('LuT_314','LuT_736'),]$disease<-'CPFE'
metadata$colnames.pseudobulk.<-NULL
rownames(metadata)<-colnames(pseudobulk)

dds<-DESeqDataSetFromMatrix(countData=pseudobulk,
                            colData=metadata,
                            design=~group)

keep<-rowSums(counts(dds)) >= 5 
dds<-dds[keep,]

rld<-rlog(dds, blind=T)                      
rld_df<-as.data.frame(assay(rld))

rv<-rowVars(assay(rld))
select<-order(rv, decreasing=TRUE)[seq_len(min(1000, length(rv)))]
pca<-prcomp(t(rld_df[select,]))
percentVar<-pca$sdev^2/sum(pca$sdev^2)
pcaData<-data.frame(PC1=pca$x[,1], PC2=pca$x[,2])
pcaData$disease<-metadata$disease[match(rownames(pcaData), rownames(metadata))]
pcaData$group<-metadata$group[match(rownames(pcaData), rownames(metadata))]
colnames(pcaData)[1:2]<-c('PC1', 'PC2')

attr(pcaData, "percentVar") <- percentVar[c(1, 2)]
percentVar<-round(100*attr(pcaData, "percentVar"), 1)

ggplot(pcaData, aes(x=PC1, y=PC2, color=group)) +
  geom_point(size=5) +
  geom_mark_ellipse(aes(fill=group), alpha=0.15) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values=c("Control"="#E6AB02",
                              "COPD"="#1B9E77",
                              "CPFE"="#D95F02",
                              "IPF"="#7570B3")) +
  scale_fill_manual(values=c("Control"="#E6AB02",
                             "COPD"="#1B9E77",
                             "CPFE"="#D95F02",
                             "IPF"="#7570B3")) +
  coord_fixed()+
  labs(title="PCA") +
  theme(panel.background=element_rect(fill=NA, colour="black"),
        aspect.ratio=0.8,
        plot.background=element_blank(),
        legend.background=element_blank(),
        plot.title=element_text(size=14, face="bold", hjust=0.5),
        axis.text=element_text(colour="black", size=12),
        axis.ticks=element_line(colour="black", size=1),
        axis.title=element_text(colour="black", size=14),
        legend.key=element_rect(fill="white", colour="black"),
        legend.text=element_text(colour="black", size=12),
        legend.title=element_text(colour="black", size=14))

# Run GSVA
list2<-list(fibrosis=c('COL1A1','COL1A2','COL3A1','COL5A1','COL5A2','COL6A1','COL6A2','COL6A3','FN1','POSTN','SPARC','TNC','MMP7',
                       'MMP3','THBS1','VCAN','COMP','LUM','DCN','BGN','COL18A1','MUC5B','DNAH6','DNAH7','DNAI1','RPGRIP1'))

es<-gsva(as.matrix(rld_df),
         list2,
         min.sz=2,
         max.sz=500,
         verbose=T,
         method = 'gsva',
         #kcdf = 'Poisson',   
         mx.diff=F,  
         parallel.sz=1)
  
es<-melt(es)
es$group<-metadata$group

ggplot(es, aes(x=group, y=value, fill=group)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("Control"="#E6AB02",
                             "COPD"="#1B9E77",
                             "CPFE"="#D95F02",
                             "IPF"="#7570B3")) +
  theme_bw() 

res.aov<-aov(value ~ group, es)
var<-leveneTest(value ~ group, es)
var<-var$`Pr(>F)`[1]
resi<-shapiro.test(residuals(object=res.aov))
resi<-resi$p.value 

if(var | resi <= 0.05){
  stat_p<-dunn.test(x=es$value, g=es$group)
  tmp<-data.frame(test='Kruskal-Wallis non parametric', comparison=stat_p$comparisons, adjusted.p=stat_p$P.adjusted)
  results<-tmp
  
  } else {
    stat_p<-TukeyHSD(res.aov)
    tmp<-data.frame(test='one-way ANOVA', comparison=rownames(stat_p$group), adjusted.p=stat_p$group[,4])
    results<-tmp
  }
}

# Load datasets
Idents(lymphnode)<-'LuT'
table(lymphnode$LuT)

# Aggregate counts
pseudobulk<-AggregateExpression(lymphnode,
                                group.by="LuT",
                                assays="RNA",
                                slot="counts")

pseudobulk<-pseudobulk$RNA
pseudobulk<-as.data.frame(pseudobulk)

metadata<-data.frame(colnames(pseudobulk))

metadata$group<-'COPD'
metadata[metadata$colnames.pseudobulk. %in% c('LuT_358','LuT_424','LuT_503','LuT_636'),]$group<-'IPF'
metadata[metadata$colnames.pseudobulk.=='LuT_737',]$group<-'CPFE'
metadata$colnames.pseudobulk.<-NULL
rownames(metadata)<-colnames(pseudobulk)

dds<-DESeqDataSetFromMatrix(countData=pseudobulk,
                            colData=metadata,
                            design=~group)

keep<-rowSums(counts(dds)) >= 5 
dds<-dds[keep,]

rld<-rlog(dds, blind=T)                      
rld_df<-as.data.frame(assay(rld))

rv<-rowVars(assay(rld))
select<-order(rv, decreasing=TRUE)[seq_len(min(1000, length(rv)))]
pca<-prcomp(t(rld_df[select,]))
percentVar<-pca$sdev^2/sum(pca$sdev^2)
pcaData<-data.frame(PC1=pca$x[,1], PC2=pca$x[,2])
pcaData$group<-metadata$group[match(rownames(pcaData), rownames(metadata))]
colnames(pcaData)[1:2]<-c('PC1', 'PC2')

attr(pcaData, "percentVar") <- percentVar[c(1, 2)]
percentVar<-round(100*attr(pcaData, "percentVar"), 1)

ggplot(pcaData, aes(x=PC1, y=PC2, color=group)) +
  geom_point(size=5) +
  geom_mark_ellipse(aes(fill=group), alpha=0.15) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values=c("Control"="#E6AB02",
                              "COPD"="#1B9E77",
                              "CPFE"="#D95F02",
                              "IPF"="#7570B3")) +
  scale_fill_manual(values=c("Control"="#E6AB02",
                             "COPD"="#1B9E77",
                             "CPFE"="#D95F02",
                             "IPF"="#7570B3")) +
  coord_fixed()+
  labs(title="PCA") +
  theme(panel.background=element_rect(fill=NA, colour="black"),
        aspect.ratio=0.8,
        plot.background=element_blank(),
        legend.background=element_blank(),
        plot.title=element_text(size=14, face="bold", hjust=0.5),
        axis.text=element_text(colour="black", size=12),
        axis.ticks=element_line(colour="black", size=1),
        axis.title=element_text(colour="black", size=14),
        legend.key=element_rect(fill="white", colour="black"),
        legend.text=element_text(colour="black", size=12),
        legend.title=element_text(colour="black", size=14))
