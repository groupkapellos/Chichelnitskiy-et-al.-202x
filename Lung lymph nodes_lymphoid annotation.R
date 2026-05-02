# Load packages
library(readxl)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# Load signatures
szabo<-list()
for (i in 1:11){
  df<-read_xlsx('Szabo 2019 Nat Commun.xlsx', skip=1, sheet=i)
  df<-df[df$fdr<=0.05, ]
  df<-df %>% arrange(-log2_effect)
  df<-df[1:100,]
  szabo[[i]]<-df$gene[!is.na(df$gene)]
}
names(szabo)<-c('CD4rest','CD4act1','CD4act2','CD4act3','CD4TRMrest','CD4TRMact','CD4Treg','CD8TRMEMrest','CD8TRMrest','CD8TRMEMact','CD8TEMRA')

naive<-c('SELL', 'PLPP5', 'FCER2', 'IL4R', 'IGHM', 'IGHD', 'TCL1A', 'MS4A1', 'CD79A', 'CD19', 'EIF4EBP1', 'TNFRSF13C',
         'CD24', 'PAX5', 'RPS27', 'VPREB1', 'CIITA', 'MYC', 'REL', 'EGFR', 'TGIF1', 'IRF8', 'TXNIP', 'FCMR', 'BANK1',
         'IL6')

mzb<-c('PLD4', 'MZB1', 'CD1C', 'SOX4', 'CXCR4', 'CD83', 'BACH2', 'CCR6', 'CD44', 'CD69', 'BANK1', ' RASGRP2', 'STAT1',
       'IFNGR1', 'CD9', 'PRDM1', 'XBP1', 'TNFRSF17', 'FKBP11')

memory<-c('TAGLN2', 'IGHA2', 'AHNAK', 'S100A4', 'CRIP2', 'ITGB1', 'VIM', 'CRIP1', 'LGALS1', 'IGHA1', 'S100A10',
          'MT-ATP8', 'CCR6', 'GPR183', 'CD44', 'KLF2', 'TNFRSF13B', 'PLAC8', 'FCRL4', 'CCR1', ' ITGAX', 'CD27',
          'COCH', 'IGHM', 'IL6', 'IL1A', 'BANK1', 'RASGRP2')

plasma<-c('IGHM', 'IGHD', 'TNFRSF13C', 'IGKC', 'IGLC2', 'MS4A1', 'CD79A', 'PRDM1', 'SDC1', 'CD27', 'XBP1', 'FOS', 
          'IRF4', 'CREB3L2', 'RASSF6', 'FRZB', 'HOPX', 'BTNL9', 'FGFR1', 'JCHAIN', 'MZB1', 'SSR4', 'TNFRSF17', 
          'IL1B', 'FKBP11')

gc<-c('CCND2', 'MIR155HG', 'PSME2', 'BHLHE40', 'PARVB', 'EBI3', 'BCL2A1', 'LMO2', 'GMDS', 'PRPSAP2', 'SERPINA9',
      'MARCKSL1', 'CD27', 'CD38', 'BCL6', 'SUGCT', 'EZR', 'ISG20', 'AICDA', 'CXCR4', 'FOXP1', 'PCNA', 'MKI67',
      'CDK1', 'CDC20', 'MME', 'CD72', 'PTPN6', 'IFNGR1', 'CAMK1', 'CD22', 'CD83', 'CD40', 'NFKBIA')

# Prepare gene sets
geneSet<-list(CD4rest=setdiff(szabo$CD4rest, unique(unlist(szabo[-1]))),
              CD4act1=setdiff(szabo$CD4act1, unique(unlist(szabo[-2]))),
              CD4act2=setdiff(szabo$CD4act2, unique(unlist(szabo[-3]))),
              CD4act3=setdiff(szabo$CD4act3, unique(unlist(szabo[-4]))),
              CD4TRMrest=setdiff(szabo$CD4TRMrest, unique(unlist(szabo[-5]))),
              CD4TRMact=setdiff(szabo$CD4TRMact, unique(unlist(szabo[-6]))),
              CD4Treg=setdiff(szabo$CD4Treg, unique(unlist(szabo[-7]))),
              CD8TRMEMrest=setdiff(szabo$CD8TRMEMrest, unique(unlist(szabo[-8]))),
              CD8TRMrest=setdiff(szabo$CD8TRMrest, unique(unlist(szabo[-9]))),
              CD8TRMEMact=setdiff(szabo$CD8TRMEMact, unique(unlist(szabo[-10]))),
              CD8TEMRA=setdiff(szabo$CD8TEMRA, unique(unlist(szabo[-11]))))

# Gene signature heatmaps
test<-tcells@assays$RNA@data  
test<-test[rowSums(as.matrix(test))>0,]

tmp2<-data.frame(rownames(test))
colnames(tmp2)<-'genes'
for (i in sort(as.character(unique(tcells@meta.data$integrated_snn_res.0.4)),decreasing=F)) {
  tmp<-test[,colnames(test) %in% rownames(tcells@meta.data[tcells@meta.data$integrated_snn_res.0.4==i,])]
  tmp<-as.data.frame(tmp)
  tmp<-as.data.frame(rowMeans(tmp))
  colnames(tmp)<-i
  tmp2<-cbind(tmp2,tmp)
}  
tmp2$genes<-NULL

test2<-tmp2[rownames(tmp2) %in% unname(unlist(geneSet)),]
test2<-test2[unname(unlist(geneSet)),]
test2<-test2[!grepl('NA', rownames(test2)),]

pheatmap(test2,
         clustering_distance_rows='correlation',
         clustering_distance_cols='correlation',
         scale='row',
         show_rownames=F,
         cluster_rows=F,
         cluster_cols=F,
         cellheight=0.5,
         cellwidth=16,
         gaps_row=c(34,103,150,167,205,242,272,316,378),
         breaks=seq(-2, 2, by=0.05),
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.05))))
