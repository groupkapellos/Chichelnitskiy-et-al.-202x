# Load packages
library(readxl)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# Gene signature heatmaps
test<-nks@assays$RNA@data  
test<-test[rowSums(as.matrix(test))>0,]

tmp2<-data.frame(rownames(test))
colnames(tmp2)<-'genes'

for (i in levels(nks$integrated_snn_res.0.3)) {
  tmp<-test[,colnames(test) %in% rownames(nks@meta.data[nks$integrated_snn_res.0.3==i,])]
  tmp<-as.data.frame(tmp)
  tmp<-as.data.frame(rowMeans(tmp))
  colnames(tmp)<-i
  tmp2<-cbind(tmp2,tmp)
}  
tmp2$genes<-NULL

features=c('TMSB4X','IL32','CCL5','CD52','ITGA1','RBPJ','RGS1','IFNGR1','FOS','CXCR4','ID2',
           'CD69','BTG1','PABPC1','RPL13','RPL13A','IL7R',
           'SELL','CCR7','IL2RA','FOXP1',
           'LTB',
           'STMN1','MKI67','LGALS1',
           'PFN1','KLRD1','PRF1','NKG7','GZMH','RUNX3','FGFBP2','LYST','ZEB2','CD47','NDUFA5','IKZF3',
           'CD99','KLRB1','KLRF1','PLAC8','SLAMF7','KLF3','S1PR5','MYO1F','CRTAM','CD7','DUSP2',
           'DTHD1','RORA','CCL3','GZMB','CCL4','GNLY','GZMA','CD44','CD74','HLA-DRA','HLA-DRB1','HOPX','LGALS3','STAT1','IFNG','BTBD9','BHLHE40','XBP1','ANXA2')

# Snyder, Galletti, Clarke, Madissoon, Szabo, Leader

features=c('MS4A1','HLA-DRB1','CD69','LTB','GPR183','CD83','BTG1','CXCR4','ACTB',          # B                         
           'IGHM',                                                                         # IgM plasma
           'IGHG1','IGHG2','IGHG3','IGHG4',                                                # IgG plasma
           'MZB1','XBP1','SSR4','JCHAIN',                                                  # plasma
           'IGHA1','IGHA2')                                                                # IgA plasma
                                                                                                     
test2<-tmp2[rownames(tmp2) %in% unique(features),]
test2<-test2[unique(features),]
test2<-test2[!grepl('NA', rownames(test2)),]

pheatmap(test2,
         clustering_distance_rows='correlation',
         clustering_distance_cols='correlation',
         scale='row',
         show_colnames=T,
         cluster_rows=F,
         cluster_cols=F,
         cellheight=8,
         cellwidth=12,
         breaks=seq(-2, 2, by=0.05),
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.05))))

DotPlot(nks, features=c('LGALS1','HMGB2','CCL4','CCL4L2','CCL3','IFNG','HLA-DPA1','HLA-DRA','HLA-DRB1','HLA-DPB1','CD74','PRF1','GZMB','FCGR3A','XCL1'), cols='RdBu', dot.min=0.1) + RotatedAxis()
