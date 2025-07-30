# Remove existing variables from R memory
rm(list=ls())

# Load packages
library(Seurat)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# All myeloid cells
DefaultAssay(object=lymphnode)<-"integrated"
Idents(lymphnode)<-lymphnode$integrated_snn_res.0.3
myeloid<-subset(lymphnode, idents=c(4, 7))

Idents(myeloid)<-'doublets'
myeloid<-subset(myeloid, idents='Singlet')

Idents(myeloid)<-'general_final_annotation'
myeloid<-subset(myeloid, idents=unique(myeloid$general_final_annotation)[c(3,5,6,7,9)])

# Find variable genes
myeloid<-FindVariableFeatures(object=myeloid, selection.method="vst", nfeatures=2000)

# Scale data
myeloid<-ScaleData(object=myeloid)

# Perform linear reduction    
myeloid<-RunPCA(object=myeloid, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=myeloid, ndims=50)

# Find k-nearest neighbours
DefaultAssay(myeloid)<-'integrated'
myeloid<-FindNeighbors(object=myeloid, reduction='pca', dims=1:15)

# Run UMAP
myeloid<-RunUMAP(object=myeloid, reduction='pca', dims=1:15)
DimPlot(object=myeloid, reduction='umap', group.by='general_final_annotation', label=F, pt.size=0.1)


## Macrophages
DefaultAssay(object=myeloid)<-"integrated"
Idents(myeloid)<-myeloid$general_final_annotation
macs<-subset(myeloid, idents=unique(myeloid$general_final_annotation)[c(1,5)])

# Find variable genes
macs<-FindVariableFeatures(object=macs, selection.method="vst", nfeatures=2000)

# Scale data
macs<-ScaleData(object=macs)

# Perform linear reduction    
macs<-RunPCA(object=macs, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=macs, ndims=50)

# Find k-nearest neighbours
DefaultAssay(macs)<-'integrated'
macs<-FindNeighbors(object=macs, reduction='pca', dims=1:20)

# Set clusters
macs<-FindClusters(object=macs, resolution=0.3, group.singletons=T)

# Run UMAP
macs<-RunUMAP(object=macs, reduction='pca', dims=1:20)
DimPlot(object=macs, reduction='umap', label=F, pt.size=0.1) 

# Find cluster markers
DefaultAssay(macs)<-'RNA'
results.macs<-all.markers(object=macs, min.pct=0.20, log=0.25)
results.macs<-results.macs[results.macs$p_val_adj<=0.05 & results.macs$avg_log2FC>0,]
results.macs<-results.macs %>% arrange(cluster, -avg_log2FC)
results.macs5<-results.macs %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)

DotPlot(object=macs, features=unique(results.macs5$gene), cols='RdBu', dot.scale=6, dot.min=0.1) + RotatedAxis()

# Plot DE genes between clusters
macs@meta.data[macs@meta.data$disease=='Fibrosis EAA',]$disease<-'Fibrosis'
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

test2<-tmp2[rownames(tmp2) %in% unique(results.macs5$gene),]
test2<-test2[unique(results.macs5$gene),]
test2<-test2[,c(1,6,
                2,7,
                3,
                4,
                5)]

df<-data.frame(disease=c(rep(c('Emphysema', 'Fibrosis'), 2),rep('Emphysema',3)),
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
         gaps_col=c(2,4,5,6),
         annotation_col=df,
         breaks=seq(-2, 2, by=0.05),
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.05))))

## Monocytes
DefaultAssay(object=myeloid)<-"integrated"
Idents(myeloid)<-myeloid$general_final_annotation
monos<-subset(myeloid, idents=unique(myeloid$general_final_annotation)[3])

# Find variable genes
monos<-FindVariableFeatures(object=monos, selection.method="vst", nfeatures=2000)

# Scale data
monos<-ScaleData(object=monos)

# Perform linear reduction    
monos<-RunPCA(object=monos, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=monos, ndims=50)

# Find k-nearest neighbours
DefaultAssay(monos)<-'integrated'
monos<-FindNeighbors(object=monos, reduction='pca', dims=1:30)

# Set clusters
monos<-FindClusters(object=monos, resolution=0.3, group.singletons=T)

# Run UMAP
monos<-RunUMAP(object=monos, reduction='pca', dims=1:30)
DimPlot(object=monos, reduction='umap', label=F, pt.size=0.3)+scale_color_igv()

# Find cluster markers
DefaultAssay(monos)<-'RNA'
results.monos<-all.markers(object=monos, min.pct=0.20, log=0.25)
results.monos<-results.monos[results.monos$p_val_adj<=0.05 & results.monos$avg_log2FC>0,]
results.monos<-results.monos %>% arrange(cluster, -avg_log2FC)
results.monos5<-results.monos %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)

DotPlot(object=monos, features=unique(results.monos5$gene), cols='RdBu', dot.scale=6, dot.min=0.1) + RotatedAxis()

# Plot DE genes between clusters
monos@meta.data[monos@meta.data$disease=='Fibrosis EAA',]$disease<-'Fibrosis'
test<-monos@assays$RNA@data  
test<-test[rowSums(as.matrix(test))>0,]

tmp2<-data.frame(rownames(test))
colnames(tmp2)<-'genes'
monos@meta.data$heatmap<-paste(monos@meta.data$disease, sep='_', monos@meta.data$integrated_snn_res.0.3)

for (i in sort(as.character(unique(monos@meta.data$heatmap)),decreasing=F)) {
  tmp<-test[,colnames(test) %in% rownames(monos@meta.data[monos@meta.data$heatmap==i,])]
  tmp<-as.data.frame(tmp)
  tmp<-as.data.frame(rowMeans(tmp))
  colnames(tmp)<-i
  tmp2<-cbind(tmp2,tmp)
}  
tmp2$genes<-NULL

test2<-tmp2[rownames(tmp2) %in% unique(results.monos5$gene),]
test2<-test2[unique(results.monos5$gene),]
test2<-test2[,c(1,4,
                2,
                3)]

df<-data.frame(disease=c('Emphysema', 'Fibrosis', rep('Emphysema', 2)),
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
         gaps_col=c(2,3),
         annotation_col=df,
         breaks=seq(-2, 2, by=0.05),
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.05))))

## Dendritic cells
DefaultAssay(object=myeloid)<-"integrated"
Idents(myeloid)<-myeloid$general_final_annotation
dcs<-subset(myeloid, idents=unique(myeloid$general_final_annotation)[c(2,4)])

# Find variable genes
dcs<-FindVariableFeatures(object=dcs, selection.method="vst", nfeatures=2000)

# Scale data
dcs<-ScaleData(object=dcs)

# Perform linear reduction    
dcs<-RunPCA(object=dcs, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=dcs, ndims=50)

# Find k-nearest neighbours
DefaultAssay(dcs)<-'integrated'
dcs<-FindNeighbors(object=dcs, reduction='pca', dims=1:25)

# Set clusters
dcs<-FindClusters(object=dcs, resolution=0.5, group.singletons=T)

# Run UMAP
dcs<-RunUMAP(object=dcs, reduction='pca', dims=1:25)
DimPlot(object=dcs, reduction='umap', label=F, pt.size=0.3)+scale_color_jama()

# Find cluster markers
DefaultAssay(dcs)<-'RNA'
results.dcs<-all.markers(object=dcs, min.pct=0.20, log=0.25)
results.dcs<-results.dcs[results.dcs$p_val_adj<=0.05 & results.dcs$avg_log2FC>0,]
results.dcs<-results.dcs %>% arrange(cluster, -avg_log2FC)
results.dcs5<-results.dcs %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)

DotPlot(object=dcs, features=unique(results.dcs5$gene), cols='RdBu', dot.scale=6, dot.min=0.1) + RotatedAxis()

# Plot DE genes between clusters
dcs@meta.data[dcs@meta.data$disease=='Fibrosis EAA',]$disease<-'Fibrosis'
test<-dcs@assays$RNA@data  
test<-test[rowSums(as.matrix(test))>0,]

tmp2<-data.frame(rownames(test))
colnames(tmp2)<-'genes'
dcs@meta.data$heatmap<-paste(dcs@meta.data$disease, sep='_', dcs@meta.data$integrated_snn_res.0.5)

for (i in sort(as.character(unique(dcs@meta.data$heatmap)),decreasing=F)) {
  tmp<-test[,colnames(test) %in% rownames(dcs@meta.data[dcs@meta.data$heatmap==i,])]
  tmp<-as.data.frame(tmp)
  tmp<-as.data.frame(rowMeans(tmp))
  colnames(tmp)<-i
  tmp2<-cbind(tmp2,tmp)
}  
tmp2$genes<-NULL

test2<-tmp2[rownames(tmp2) %in% unique(results.dcs5$gene),]
test2<-test2[unique(results.dcs5$gene),]
test2<-test2[,c(1,3,
                2,4,
                5)]

df<-data.frame(disease=c(rep(c('Emphysema', 'Fibrosis'), 2), 'Fibrosis'),
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
         gaps_col=c(2,4),
         annotation_col=df,
         breaks=seq(-2, 2, by=0.05),
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.05))))