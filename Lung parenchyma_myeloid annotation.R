# Load packages
library(readxl)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# Load signatures
morse<-read_xlsx('Morse 2019 ERJ.xlsx')
morse<-morse[morse$p_adjust<=0.05,]
morse<-morse[morse$Logfoldchange>=1,]  
morse<-morse[morse$Cluster %in% c(0,1,3),]

day36<-read_xlsx('Evren 2021 Immunity.xlsx', skip=2)
day36<-day36[day36$p_val_adj<=0.05,]
day36<-day36[day36$cluster %in% c(0,1,2,3,5),]

florent<-read_xlsx('Dutertre 2019 Immunity.xlsx', sheet=1)
florent$p_val_adj<-as.numeric(florent$p_val_adj)
florent<-florent[florent$p_val_adj<=0.05,]
florent<-florent[!is.na(florent$gene),]
florent<-florent[florent$cluster %in% c(2,4,6,8),]
florent$gene<-sapply(strsplit(split='.', florent$gene, fixed=T), '[', 1)
florent<-florent %>% arrange(cluster, -avg_logFC)

# Prepare gene sets
mor<-list()
for (i in unique(morse$Cluster)){
  tmp<-setdiff(morse[morse$Cluster==i,]$Gene, unique(morse[morse$Cluster!=i,]$Gene))
  mor[[i]]<-tmp
}
names(mor)<-c('FABP4 macs', 'SPP1 macs', 'FN1 monos')

flo<-list()
for (i in c(2,4,6,8)){
  flo[[i]]<-florent[florent$cluster==i,]$gene
}
flo<-flo[c(2,4,6,8)]

flo2<-list(setdiff(flo[[1]], unique(c(flo[[2]],flo[[3]],flo[[4]]))),
           setdiff(flo[[2]], unique(c(flo[[1]],flo[[3]],flo[[4]]))),
           setdiff(flo[[3]], unique(c(flo[[2]],flo[[1]],flo[[4]]))),
           setdiff(flo[[4]], unique(c(flo[[2]],flo[[3]],flo[[1]]))))

names(flo2)<-c('cDC2/pre', 'cDC2', 'cDC1', 'pDC')

# Gene signature heatmaps
test<-dcs@assays$RNA@data  
test<-test[rowSums(as.matrix(test))>0,]

tmp2<-data.frame(rownames(test))
colnames(tmp2)<-'genes'

for (i in sort(as.character(unique(dcs@meta.data$integrated_snn_res.0.2)),decreasing=F)) {
  tmp<-test[,colnames(test) %in% rownames(dcs@meta.data[dcs@meta.data$integrated_snn_res.0.2==i,])]
  tmp<-as.data.frame(tmp)
  tmp<-as.data.frame(rowMeans(tmp))
  colnames(tmp)<-i
  tmp2<-cbind(tmp2,tmp)
}  
tmp2$genes<-NULL

genes<-c(mor[[1]], mor[[2]], mor[[3]])
genes<-c('MS4A6A','LGALS2','CPVL','CD69','HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQB1',
         'CTNNB1','C1orf56','CDC42SE1',
         'FCGR3A','LST1','AIF1','CD52',
         'S100A8','S100A9','S100A12','VCAN','LYZ','CD14',
         'MALAT1','NEAT1','MT-CYB','MT-CO1','MT-CO2','MT-CO3','MT-ND2','MT-ND3','MT-ND4','MT-ND5','MT-ATP6','UCP2')
genes<-c(flo2[[1]],flo2[[2]],flo2[[3]],flo2[[4]])

test2<-tmp2[rownames(tmp2) %in% genes,]
test2<-test2[genes,]
test2<-test2[!grepl('NA', rownames(test2)),]

pheatmap(test2,
         clustering_distance_rows='correlation',
         clustering_distance_cols='correlation',
         scale='row',
         show_colnames=T,
         show_rownames=F,
         cluster_rows=F,
         cluster_cols=F,
         cellheight=3,
         cellwidth=16,
         gaps_row=c(13,14,74),
         breaks=seq(-2, 2, by=0.05),
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.05))))
