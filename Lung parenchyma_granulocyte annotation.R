# Remove unneeded objects 
rm(list=ls())

# Load packages
library(readxl)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# Load signatures
zilionis<-read_xlsx('Zilionis 2019 Immunity.xlsx', sheet=4)
zilionis<-as.data.frame(zilionis)

zilio<-list()
for (i in 2:6){
  tmp<-zilionis[order(zilionis[,i], decreasing=T),]
  zilio[[i]]<-tmp[tmp[,i]>=0.99,]$...1
}
zilio<-zilio[2:6]
names(zilio)<-colnames(zilionis)[2:6]

trava<-list()
for (i in 1:2){
  tmp<-read_xlsx('Travaglini 2020 Nature.xlsx', sheet=i, skip=1)
  tmp<-as.data.frame(tmp)
  tmp<-tmp[tmp$avg_logFC>=1 & tmp$p_val_adj<=0.05,]
  trava[[i]]<-tmp$Gene
}

# Prepare gene sets
zilio2<-list(setdiff(zilio[[1]], unique(c(zilio[[2]],zilio[[3]],zilio[[4]],zilio[[5]]))),
             setdiff(zilio[[2]], unique(c(zilio[[1]],zilio[[3]],zilio[[4]],zilio[[5]]))),
             setdiff(zilio[[3]], unique(c(zilio[[2]],zilio[[1]],zilio[[4]],zilio[[5]]))),
             setdiff(zilio[[4]], unique(c(zilio[[2]],zilio[[3]],zilio[[1]],zilio[[5]]))),
             setdiff(zilio[[5]], unique(c(zilio[[2]],zilio[[3]],zilio[[4]],zilio[[1]]))))
             
names(zilio2)<-colnames(zilionis)[2:6]

trava2<-list(setdiff(trava[[1]], trava[[2]]),
             setdiff(trava[[2]], trava[[1]]))

# Gene signature heatmaps
test<-neutros@assays$RNA@data  
test<-test[rowSums(as.matrix(test))>0,]

tmp2<-data.frame(rownames(test))
colnames(tmp2)<-'genes'

for (i in sort(as.character(unique(neutros@meta.data$integrated_snn_res.0.4)),decreasing=F)) {
  tmp<-test[,colnames(test) %in% rownames(neutros@meta.data[neutros@meta.data$integrated_snn_res.0.4==i,])]
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