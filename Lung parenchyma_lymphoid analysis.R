# Remove existing variables from R memory
rm(list=ls())

# Load packages
library(Seurat)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# All lymphoid cells
DefaultAssay(object=combined_list)<-"integrated"
Idents(combined_list)<-combined_list$integrated_snn_res.0.4
lymphoid<-subset(combined_list, idents=c(0,2,9,13,14,15))

Idents(lymphoid)<-'doublets_individual'
lymphoid<-subset(lymphoid, idents='Singlet')

Idents(lymphoid)<-'general_final_annotation'
lymphoid<-subset(lymphoid, idents=unique(lymphoid$general_final_annotation)[c(3:9)])

# Find variable genes
lymphoid<-FindVariableFeatures(object=lymphoid, selection.method="vst", nfeatures=2000)

# Scale data
lymphoid<-ScaleData(object=lymphoid)

# Perform linear reduction    
lymphoid<-RunPCA(object=lymphoid, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=lymphoid, ndims=50)

# Find k-nearest neighbours
DefaultAssay(lymphoid)<-'integrated'
lymphoid<-FindNeighbors(object=lymphoid, reduction='pca', dims=1:15)

# Run UMAP
lymphoid<-RunUMAP(object=lymphoid, reduction='pca', dims=1:15)
DimPlot(object=lymphoid, reduction='umap', group.by='general_final_annotation', label=F, pt.size=0.1)+scale_color_manual(labels=unique(lymphoid$general_final_annotation), values=c('#309343','#B89B74','#865FAB','#D4A55B','#ED4F50','#B294C7','#E4201F'))


## T cells
DefaultAssay(object=lymphoid)<-"integrated"
Idents(lymphoid)<-lymphoid$general_final_annotation
tcells<-subset(lymphoid, idents=c('T cells 1','T cells 2','Cycling T/NK cells'))

# Find variable genes
tcells<-FindVariableFeatures(object=tcells, selection.method="vst", nfeatures=2000)

# Scale data
tcells<-ScaleData(object=tcells)

# Perform linear reduction    
tcells<-RunPCA(object=tcells, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=tcells, ndims=50)

# Find k-nearest neighbours
DefaultAssay(tcells)<-'integrated'
tcells<-FindNeighbors(object=tcells, reduction='pca', dims=1:20)

# Find clusters
tcells<-FindClusters(object=tcells, resolution=0.7, group.singletons=T)

# Run UMAP
tcells<-RunUMAP(object=tcells, reduction='pca', dims=1:20)
DimPlot(object=tcells, reduction='umap', label=F, pt.size=0.1)+scale_color_d3(palette='category20')

# Find cluster markers
DefaultAssay(tcells)<-'RNA'
results.tcells<-all.markers(object=tcells, min.pct=0.20, log=0.25)
results.tcells<-results.tcells[results.tcells$p_val_adj<=0.05 & results.tcells$avg_log2FC>0,]
results.tcells<-results.tcells %>% arrange(cluster, -avg_log2FC)
results.tcells5<-results.tcells %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)

DotPlot(object=tcells, features=unique(results.tcells5$gene), cols='RdBu', dot.scale=6, dot.min=0.1) + RotatedAxis()

# Plot DE genes between clusters
tcells@meta.data[tcells@meta.data$disease=='Fibrosis EAA',]$disease<-'Fibrosis'
test<-tcells@assays$RNA@data  
test<-test[rowSums(as.matrix(test))>0,]

tmp2<-data.frame(rownames(test))
colnames(tmp2)<-'genes'
tcells@meta.data$heatmap<-paste(tcells@meta.data$disease, sep='_', tcells@meta.data$integrated_snn_res.0.7)

for (i in sort(as.character(unique(tcells@meta.data$heatmap)),decreasing=F)) {
  tmp<-test[,colnames(test) %in% rownames(tcells@meta.data[tcells@meta.data$heatmap==i,])]
  tmp<-as.data.frame(tmp)
  tmp<-as.data.frame(rowMeans(tmp))
  colnames(tmp)<-i
  tmp2<-cbind(tmp2,tmp)
}  
tmp2$genes<-NULL

test2<-tmp2[rownames(tmp2) %in% unique(results.tcells5$gene),]
test2<-test2[unique(results.tcells5$gene),]
test2<-test2[,c(28,1,11,20,
                29,2,12,21,
                31,3,13,22,
                32,4,14,23,
                33,5,24,
                34,6,16,25,
                35,7,17,26,
                36,8,18,27,
                37,9,
                10,
                30)]

df<-data.frame(disease=c(rep(c('Tumor-free', 'Emphysema', 'Fibrosis', 'PAH'), 4),'Tumor-free', 'Emphysema', 'PAH',rep(c('Tumor-free', 'Emphysema', 'Fibrosis', 'PAH'), 3),'Tumor-free', 'Emphysema', 'Emphysema','Tumor-free'),
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
         gaps_col=c(4,8,12,16,19,23,27,31,33,34),
         annotation_col=df,
         breaks=seq(-2, 2, by=0.05),
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.05))))

## B cells
DefaultAssay(object=lymphoid)<-"integrated"
Idents(lymphoid)<-lymphoid$general_final_annotation
bcells<-subset(lymphoid, idents=c('B cells','Plasma cells'))

# Find variable genes
bcells<-FindVariableFeatures(object=bcells, selection.method="vst", nfeatures=2000)

# Scale data
bcells<-ScaleData(object=bcells)

# Perform linear reduction    
bcells<-RunPCA(object=bcells, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=bcells, ndims=50)

# Find k-nearest neighbours
DefaultAssay(bcells)<-'integrated'
bcells<-FindNeighbors(object=bcells, reduction='pca', dims=1:16)

# Find clusters
bcells<-FindClusters(object=bcells, resolution=0.2, group.singletons=T)

# Run UMAP
bcells<-RunUMAP(object=bcells, reduction='pca', dims=1:16)
DimPlot(object=bcells, reduction='umap', label=F, pt.size=0.1)+scale_color_ucscgb()

# Find cluster markers
DefaultAssay(bcells)<-'RNA'
results.bcells<-all.markers(object=bcells, min.pct=0.20, log=0.25)
results.bcells<-resultsr.bcells[results.bcells$p_val_adj<=0.05 & results.bcells$avg_log2FC>0,]
results.bcells<-results.bcells %>% arrange(cluster, -avg_log2FC)
results.bcells5<-results.bcells %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)

DotPlot(object=bcells, features=unique(results.bcells5$gene), cols='RdBu', dot.scale=6, dot.min=0.1) + RotatedAxis()

# Plot DE genes between clusters
bcells@meta.data[bcells@meta.data$disease=='Fibrosis EAA',]$disease<-'Fibrosis'
test<-bcells@assays$RNA@data  
test<-test[rowSums(as.matrix(test))>0,]

tmp2<-data.frame(rownames(test))
colnames(tmp2)<-'genes'
bcells@meta.data$heatmap<-paste(bcells@meta.data$disease, sep='_', bcells@meta.data$integrated_snn_res.0.2)

for (i in sort(as.character(unique(bcells@meta.data$heatmap)),decreasing=F)) {
  tmp<-test[,colnames(test) %in% rownames(bcells@meta.data[bcells@meta.data$heatmap==i,])]
  tmp<-as.data.frame(tmp)
  tmp<-as.data.frame(rowMeans(tmp))
  colnames(tmp)<-i
  tmp2<-cbind(tmp2,tmp)
}  
tmp2$genes<-NULL

test2<-tmp2[rownames(tmp2) %in% unique(results.bcells5$gene),]
test2<-test2[unique(results.bcells5$gene),]
test2<-test2[,c(14,1,5,10,
                15,2,6,11,
                16,3,7,12,
                17,4,8,13,
                9)]

df<-data.frame(disease=c(rep(c('Tumor-free', 'Emphysema', 'Fibrosis', 'PAH'), 4),'Emphysema'),
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
         gaps_col=c(4,8,12,16),
         annotation_col=df,
         breaks=seq(-2, 2, by=0.05),
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.05))))

## NK cells
DefaultAssay(object=lymphoid)<-"integrated"
Idents(lymphoid)<-lymphoid$general_final_annotation
nks<-subset(lymphoid, idents=c('NK cells', 'CD56 NK cells'))

# Find variable genes
nks<-FindVariableFeatures(object=nks, selection.method="vst", nfeatures=2000)

# Scale data
nks<-ScaleData(object=nks)

# Perform linear reduction    
nks<-RunPCA(object=nks, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=nks, ndims=50)

# Find k-nearest neighbours
DefaultAssay(nks)<-'integrated'
nks<-FindNeighbors(object=nks, reduction='pca', dims=1:20)

# Find clusters
nks<-FindClusters(object=nks, resolution=0.3, group.singletons=T)

# Run UMAP
nks<-RunUMAP(object=nks, reduction='pca', dims=1:20)
DimPlot(object=nks, reduction='umap', label=F, pt.size=0.1)+scale_color_aaas()

# Find cluster markers
DefaultAssay(nks)<-'RNA'
results.nks<-all.markers(object=nks, min.pct=0.20, log=0.25)
results.nks<-resultsr.nks[results.nks$p_val_adj<=0.05 & results.nks$avg_log2FC>0,]
results.nks<-results.nks %>% arrange(cluster, -avg_log2FC)
results.nks5<-results.nks %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)

DotPlot(object=nks, features=unique(results.nks5$gene), cols='RdBu', dot.scale=6, dot.min=0.1) + RotatedAxis

# Plot DE genes between clusters
nks@meta.data[nks@meta.data$disease=='Fibrosis EAA',]$disease<-'Fibrosis'
test<-nks@assays$RNA@data  
test<-test[rowSums(as.matrix(test))>0,]

tmp2<-data.frame(rownames(test))
colnames(tmp2)<-'genes'
nks@meta.data$heatmap<-paste(nks@meta.data$disease, sep='_', nks@meta.data$integrated_snn_res.0.3)

for (i in sort(as.character(unique(nks@meta.data$heatmap)),decreasing=F)) {
  tmp<-test[,colnames(test) %in% rownames(nks@meta.data[nks@meta.data$heatmap==i,])]
  tmp<-as.data.frame(tmp)
  tmp<-as.data.frame(rowMeans(tmp))
  colnames(tmp)<-i
  tmp2<-cbind(tmp2,tmp)
}  
tmp2$genes<-NULL

test2<-tmp2[rownames(tmp2) %in% unique(results.nks5$gene),]
test2<-test2[unique(results.nks5$gene),]
test2<-test2[,c(13,1,5,9,
                14,2,6,10,
                15,3,7,11,
                16,4,8,12)]

df<-data.frame(disease=rep(c('Tumor-free', 'Emphysema', 'Fibrosis', 'PAH'), 4),
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
         gaps_col=c(4,8,12,16),
         annotation_col=df,
         breaks=seq(-2, 2, by=0.05),
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.05))))