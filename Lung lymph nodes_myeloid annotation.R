# Load packages
library(readxl)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# Load signatures
kim<-read_xlsx('Kim 2020 Nat Comm.xlsx')
kim<-as.data.frame(kim)

huang<-read_xlsx('Huang 2021 Cell.xlsx', sheet=1, skip=3)
huang<-huang %>% arrange(`Cell type`, -`Average log fold change`)
huang<-as.data.frame(huang)
huang<-huang[huang$`Cell type` %in% c("cDC1 (MacDC0)", "Migratory DCs (MacDC1)", "cDC2 (MacDC2)", "Macrophages (MacDC3)", "Aire APCs (MacDC4)", "pDC"),]

mouse_human_genes<-read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
mouse<-split.data.frame(mouse_human_genes, mouse_human_genes$Common.Organism.Name)[[2]]
human<-split.data.frame(mouse_human_genes, mouse_human_genes$Common.Organism.Name)[[1]]
mouse<-mouse[,c(1,4)]
human<-human[,c(1,4)]
mh_data<-merge.data.frame(mouse, human, by="DB.Class.Key", all.y=TRUE) 

huang$human<-mh_data$Symbol.y[match(huang$Gene, mh_data$Symbol.x)]

# Prepare gene sets
ki<-list()
for (i in unique(kim$Cluster)){
  tmp<-kim[kim$Cluster==i,]$Genes[1:50]
  ki[[i]]<-tmp[!is.na(tmp)]
}

hua<-list()
for (i in unique(huang$`Cell type`)[c(1,3,4,5,6)]){
  tmp<-huang[huang$`Cell type`==i,]$human[1:50]
  hua[[i]]<-tmp[!is.na(tmp)]
}

# Gene signature heatmaps
test<-monos@assays$RNA@data  
test<-test[rowSums(as.matrix(test))>0,]

tmp2<-data.frame(rownames(test))
colnames(tmp2)<-'genes'

for (i in sort(as.character(unique(monos@meta.data$integrated_snn_res.0.3)),decreasing=F)) {
  tmp<-test[,colnames(test) %in% rownames(monos@meta.data[monos@meta.data$integrated_snn_res.0.3==i,])]
  tmp<-as.data.frame(tmp)
  tmp<-as.data.frame(rowMeans(tmp))
  colnames(tmp)<-i
  tmp2<-cbind(tmp2,tmp)
}  
tmp2$genes<-NULL

genes<-c(setdiff(ki$`LN-Mac-S1`, unique(unlist(ki[-1]))),
         setdiff(ki$`LN-Mac-S2`, unique(unlist(ki[-2]))),
         setdiff(ki$`LN-Mac-S3`, unique(unlist(ki[-3]))),
         setdiff(ki$`LN-Mac-S4`, unique(unlist(ki[-4]))),
         setdiff(ki$`LN-Mac-S5`, unique(unlist(ki[-5]))),
         setdiff(ki$`LN-Mac-S6`, unique(unlist(ki[-6]))),
         setdiff(ki$`LN-Mac-S7`, unique(unlist(ki[-7]))))

genes<-c(setdiff(hua$`Aire APCs (MacDC4)`, unique(unlist(hua[-1]))),
         setdiff(hua$`cDC1 (MacDC0)`, unique(unlist(hua[-3]))),
         setdiff(hua$`cDC2 (MacDC2)`, unique(unlist(hua[-4]))),
         setdiff(hua$`Migratory DCs (MacDC1)`, unique(unlist(hua[-2]))),
         setdiff(hua$pDC, unique(unlist(hua[-5]))))

test2<-tmp2[rownames(tmp2) %in% genes,]
test2<-test2[genes,]
test2<-test2[!grepl('NA', rownames(test2)),]

pheatmap(test2,
         clustering_distance_rows='correlation',
         clustering_distance_cols='correlation',
         scale='row',
         show_colnames=T,
         cluster_rows=F,
         cluster_cols=F,
         cellheight=2,
         cellwidth=16,
         gaps_row=c(41,42,52,66,80,113),
         breaks=seq(-2, 2, by=0.05),
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.05))))
