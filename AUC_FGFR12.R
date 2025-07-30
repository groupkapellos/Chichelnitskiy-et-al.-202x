# Remove unneeded objects 
rm(list=ls())

# Set directory
setwd('~/hannover/analysis/Theo_analysis')

# Load packages
library(Seurat)
library(readxl)
library(dplyr)
library(biomaRt)
library(OmnipathR)
library(AUCell)
library(reshape2)
library(doMC)
library(doRNG)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)

# Plot FGFR expression in immune cells
dot<-DotPlot(macs, features=c('FGF2','FGFR1','FGFR2','FGFR3'), group.by='integrated_snn_res.0.3', split.by='disease', cols='RdBu') + RotatedAxis()
dot<-DotPlot(macs, features=c('FGF2','FGFR1','FGFR2','FGFR3'), group.by='integrated_snn_res.0.4', split.by='disease', cols='RdBu') + RotatedAxis()
dot<-DotPlot(macs, features=c('FGF2','FGFR1','FGFR2','FGFR3'), group.by='integrated_snn_res.0.2', split.by='disease', cols='RdBu') + RotatedAxis()
dot<-DotPlot(macs, features=c('FGF2','FGFR1','FGFR2','FGFR3'), group.by='integrated_snn_res.0.4', split.by='disease', cols='RdBu') + RotatedAxis()
dot<-DotPlot(mast, features=c('FGF2','FGFR1','FGFR2','FGFR3'), group.by='integrated_snn_res.0.2', split.by='disease', cols='RdBu') + RotatedAxis()
dot<-DotPlot(macs, features=c('FGF2','FGFR1','FGFR2','FGFR3'), group.by='integrated_snn_res.0.7', split.by='disease', cols='RdBu') + RotatedAxis()
dot<-DotPlot(macs, features=c('FGF2','FGFR1','FGFR2','FGFR3'), group.by='integrated_snn_res.0.2', split.by='disease', cols='RdBu') + RotatedAxis()
dot<-DotPlot(macs, features=c('FGF2','FGFR1','FGFR2','FGFR3'), group.by='integrated_snn_res.0.3', split.by='disease', cols='RdBu') + RotatedAxis()

dot<-DotPlot(macs, features=c('FGF2','FGFR1','FGFR2'), group.by='integrated_snn_res.0.3', split.by='disease', cols='RdBu') + RotatedAxis()
dot<-DotPlot(macs, features=c('FGF2','FGFR1','FGFR2'), group.by='integrated_snn_res.0.3', split.by='disease', cols='RdBu') + RotatedAxis()
dot<-DotPlot(macs, features=c('FGF2','FGFR1','FGFR2'), group.by='integrated_snn_res.0.5', split.by='disease', cols='RdBu') + RotatedAxis()
dot<-DotPlot(macs, features=c('FGF2','FGFR1','FGFR2'), group.by='integrated_snn_res.0.4', split.by='disease', cols='RdBu') + RotatedAxis()
dot<-DotPlot(macs, features=c('FGF2','FGFR1','FGFR2'), group.by='integrated_snn_res.0.2', split.by='disease', cols='RdBu') + RotatedAxis()
dot<-DotPlot(macs, features=c('FGF2','FGFR1','FGFR2'), group.by='integrated_snn_res.0.3', split.by='disease', cols='RdBu') + RotatedAxis()

# Assess overexpression in fibrosis
dot<-as.data.frame(dot$data)
dot[]<-lapply(dot, function(x) if(is.factor(x)) as.character(x) else x)
dot$cluster<-sapply(strsplit(split='_', dot$id), '[', 1)
dot$disease<-sapply(strsplit(split='_', dot$id), '[', 2)

dot<-dot[dot$cluster %in% c(0,1,2,3,4,5,6),]
dot<-dot[dot$disease %in% c('Emphysema', 'Fibrosis'),]

result_list<-list()
for (i in unique(dot$cluster)) {
  
  df<-dot[dot$cluster==i,]
  
  tmp2<-data.frame()
  
  for (j in unique(df$features.plot)) {
    
      avg_exp_emphysema<-df[df$features.plot == j & df$disease == 'Emphysema', "avg.exp"]
      avg_exp_fibrosis<-df[df$features.plot == j & df$disease == 'Fibrosis', "avg.exp"]
      
      if (length(avg_exp_emphysema) > 0 && length(avg_exp_fibrosis) > 0 && avg_exp_emphysema > 0) {
        if (avg_exp_fibrosis / avg_exp_emphysema > 1) {
          tmp<-data.frame(cluster=i, ligand=j)
          tmp2<-rbind(tmp2, tmp)        
        }
      }
    }
  result_list[[i]]<-tmp2
}
dot<-do.call(rbind, result_list)

# Rank genes
tmp<-as.matrix(tcells@assays$RNA@data)
tmp2<-rowSums(tmp)[rowSums(tmp)>3]
tmp<-tmp[row.names(tmp) %in% names(tmp2),]
dim(tmp)

cells_rankings<-AUCell_buildRankings(tmp)
cells_rankings

# Prepare gene sets
downstream<-import_pathwayextra_interactions()
downstream<-downstream[downstream$source_genesymbol=='FGFR2',]
downstream<-unique(downstream$target_genesymbol)

# Run AUCell analysis
cells_AUC<-AUCell_calcAUC(downstream, cells_rankings, aucMaxRank=ceiling(0.05*nrow(cells_rankings)), normAUC=T)
cells_AUC

# Plot result
res<-as.data.frame(t(getAUC(cells_AUC)))
res$cluster<-tcells$integrated_snn_res.0.4
res$disease<-tcells$disease
res<-melt(res, id.vars=c('cluster','disease'))
res<-res[res$disease %in% c('Emphysema', 'Fibrosis'),]
res<-res[res$cluster %in% c(0),]

res %>% group_by(cluster, disease) %>% summarize(Mean=mean(value))->df.tmp
df.tmp$Mean<-scale(df.tmp$Mean)
res<-merge(res, df.tmp, key="disease", all=T)

# Perform statistical analysis
#df<-data.frame()
  
for(i in unique(as.character(res$cluster))){
  tmp<-res[res$cluster==i,]
  
  a<-shapiro.test(tmp[tmp$disease=='Emphysema',]$value)
  b<-shapiro.test(tmp[tmp$disease=='Fibrosis',]$value)
  var<-var.test(value ~ disease, data=tmp)
  
  if(tmp[tmp$disease=='Fibrosis',]$Mean[1] > tmp[tmp$disease=='Emphysema',]$Mean[1]){
    print(i)
    
    if(a$p.value | b$p.value <= 0.05) {
      result<-wilcox.test(tmp[tmp$disease=='Emphysema',]$value, tmp[tmp$disease=='Fibrosis',]$value, alternative="two.sided")
      print(result$p.value)
    } else if(var$p.value <= 0.05) {
      result<-t.test(value ~ disease, data=tmp, var.equal=FALSE)
      print(result$p.value)
    } else {
      result<-t.test(value ~ disease, data=tmp, var.equal=TRUE)
      print(result$p.value)
    }
    df2<-data.frame(celltype='nks', cluster=i, gene='FGFR2', pvalue=result$p.value)
    df<-rbind(df, df2)
  } else {
    print('Signature less in fibrosis')
  }
}

#df<-rbind(df, c('neutros', 0, 'FGFR1', 1))
#df<-rbind(df, c('tcells', 1, 'FGFR1', 1))
#df<-rbind(df, c('tcells', 3, 'FGFR1', 1))
#df<-rbind(df, c('dcs', 2, 'FGFR2', 1))
#df<-rbind(df, c('neutros', 1, 'FGFR2', 1))
#df<-rbind(df, c('tcells', 5, 'FGFR2', 1))

#df<-rbind(df, c('bcells', 0, 'FGFR1', 1))
#df<-rbind(df, c('bcells', 1, 'FGFR1', 1))

df$merge<-paste(df$celltype, sep='_', df$cluster)
df$pvalue<-as.numeric(df$pvalue)

ggplot(df, aes(x=merge ,y=gene, fill=pvalue)) +
  geom_tile() +
  scale_fill_distiller(palette='Reds', breaks=c(0,0.05,0.25,0.5,0.75,1)) +
  theme(text=element_text(size=14),
        axis.title.y=element_text(size=14, face='bold'),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=14, face='bold'),
        axis.text.y=element_text(size=14, colour='black'),
        strip.text=element_text(size=12),
        strip.background=element_rect(fill='gray84', colour='black'),
        legend.title=element_text(size=14, face='bold'),
        legend.position='right',
        legend.key=element_rect(colour='black', fill='white'),
        legend.spacing.x=unit(0.3, 'cm'),
        panel.background=element_rect(fill='white'),
        panel.border=element_blank(),
        axis.line=element_line(colour="black"))

# Plot signature genes
test<-macs@assays$RNA@data  
test<-test[rowSums(as.matrix(test))>0,]

tmp2<-data.frame(rownames(test))
colnames(tmp2)<-'genes'
macs@meta.data$heatmap<-paste(macs@meta.data$disease, sep='_', macs@meta.data$integrated_snn_res.0.3)

for (i in sort(as.character(unique(macs@meta.data$heatmap)),decreasing=F)) {
  tmp<-test[,colnames(test) %in% rownames(macs@meta.data[macs@meta.data$heatmap==i,])]
  tmp<-as.data.frame(tmp)
  tmp<-as.data.frame(rowMeans(tmp))
  colnames(tmp)<-i
  tmp2<-cbind(tmp2,tmp)
}  
tmp2$genes<-NULL

downstream<-import_pathwayextra_interactions()
downstream<-downstream[downstream$source_genesymbol=='FGFR2',]
downstream<-unique(downstream$target_genesymbol)

test2<-tmp2[rownames(tmp2) %in% downstream,]
test2<-test2[,grepl('Emphysema|Fibrosis', colnames(test2))]
test2<-test2[,c(2,5)]

#list1<-list()

for(i in unique(sapply(strsplit(split='_', colnames(test2)), '[', 2))){
  df2<-test2[,grepl(i, colnames(test2))]
  df2$ratio<-df2[,grepl('Fibrosis', colnames(df2))] / df2[,grepl('Emphysema', colnames(df2))]
  df2<-df2[order(df2$ratio, decreasing=T),]
  df2$ranking<-seq(1, nrow(df2))
  df2$color<-ifelse(df2$ratio>=1, 'purple', 'green')
  df2$genes<-rownames(df2)
  list1[[paste0('macs_', i)]]<-df2
}

list1<-lapply(list1, function(x){x[order(x$genes),]})

columns_to_include<-'ranking'
merged_df<-list1[[1]][,c('genes', columns_to_include)]

for(i in 2:length(list1)){
  merged_df<-merge(merged_df, list1[[i]][,c('genes', columns_to_include)], by='genes', all=TRUE)
}

merged_df$sum<-rowSums(merged_df[,2:5])
merged_df<-merged_df[order(merged_df$sum, decreasing=FALSE),]

columns_to_include<-'color'
merged_df2<-list1[[1]][,c('genes', columns_to_include)]

for(i in 2:length(list1)){
  merged_df2<-merge(merged_df2, list1[[i]][,c('genes', columns_to_include)], by='genes', all=TRUE)
}
colnames(merged_df2)[2:5]<-names(list1)

merged_df2<-melt(merged_df2, id.vars='genes')
merged_df2$genes<-factor(merged_df2$genes, levels=rev(merged_df$genes))

ggplot(merged_df2, aes(x=variable ,y=genes, fill=value)) +
  geom_tile() +
  scale_fill_manual(values=c('#1B9E77','#7570B3')) +
  theme(text=element_text(size=14),
        axis.title.y=element_text(size=14, face='bold'),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=14, face='bold'),
        axis.text.y=element_text(size=14, colour='black'),
        strip.text=element_text(size=12),
        strip.background=element_rect(fill='gray84', colour='black'),
        legend.title=element_text(size=14, face='bold'),
        legend.position='right',
        legend.key=element_rect(colour='black', fill='white'),
        legend.spacing.x=unit(0.3, 'cm'),
        panel.background=element_rect(fill='white'),
        panel.border=element_blank(),
        axis.line=element_line(colour="black"))

########
df<-data.frame(disease=c(rep(c('Tumor-free', 'Emphysema', 'Fibrosis', 'PAH'), 5),'Fibrosis'),
               clusters=sapply(strsplit(split='_', colnames(test2)), '[', 2))
rownames(df)<-colnames(test2)

pheatmap(test2,
         clustering_distance_rows='correlation',
         clustering_distance_cols='correlation',
         scale='row',
         show_colnames=F,
         cluster_rows=F,
         cluster_cols=F,
         cellheight=9,
         cellwidth=12,
         gaps_col=c(4,8,12,16,20),
         annotation_col=df,
         breaks=seq(-2, 2, by=0.05),
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.05))))

