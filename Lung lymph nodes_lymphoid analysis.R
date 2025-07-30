# Remove existing variables from R memory
rm(list=ls())

# Load packages
library(Seurat)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# All lymphoid cells
DefaultAssay(object=lymphnode)<-"integrated"
Idents(lymphnode)<-lymphnode$integrated_snn_res.0.3
lymphoid<-subset(lymphnode, idents=c(0,3,1,2,9,6,5,7))

Idents(lymphoid)<-'doublets'
lymphoid<-subset(lymphoid, idents='Singlet')

Idents(lymphoid)<-'general_final_annotation'
lymphoid<-subset(lymphoid, idents=unique(lymphoid$general_final_annotation)[c(1,5,4,3,12,6,2,7,9,10)])

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
lymphoid<-FindNeighbors(object=lymphoid, reduction='pca', dims=1:13)

# Run UMAP
lymphoid<-RunUMAP(object=lymphoid, reduction='pca', dims=1:13)


## T cells
DefaultAssay(object=lymphoid)<-"integrated"
Idents(lymphoid)<-lymphoid$general_final_annotation
tcells<-subset(lymphoid, idents=c('T cells 1','T cells 2','Proliferating T cells','Proliferating T/B cells'))

Idents(tcells)<-'doublets'
tcells<-subset(tcells, idents='Singlet')

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
tcells<-FindNeighbors(object=tcells, reduction='pca', dims=1:10)

# Find clusters
tcells<-FindClusters(object=tcells, resolution=0.4, group.singletons=T)

# Run UMAP
tcells<-RunUMAP(object=tcells, reduction='pca', dims=1:10)
DimPlot(object=tcells, reduction='umap', label=F, pt.size=0.1)+scale_color_d3(palette='category20')

# Find cluster markers
DefaultAssay(object=tcells)<-"RNA"
results.tcells<-all.markers(object=tcells, min.pct=0.20, log=0.25)
results.tcells<-results.tcells[results.tcells$p_val_adj<=0.05 & results.tcells$avg_log2FC>0,]
results.tcells<-results.tcells %>% arrange(cluster, -avg_log2FC)
results.tcells5<-results.tcells %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)

DotPlot(object=tcells, features=unique(results.tcells5$gene), cols='RdBu', dot.scale=6, dot.min=0.1) + RotatedAxis()

# Plot DE genes between clusters
test<-tcells@assays$RNA@data  
test<-test[rowSums(as.matrix(test))>0,]

tmp2<-data.frame(rownames(test))
colnames(tmp2)<-'genes'
tcells@meta.data$heatmap<-paste(tcells@meta.data$disease, sep='.', tcells@meta.data$integrated_snn_res.0.4)

for (i in sort(as.character(unique(tcells@meta.data$heatmap)),decreasing=F)) {
  tmp<-test[,colnames(test) %in% rownames(tcells@meta.data[tcells@meta.data$heatmap==i,])]
  tmp<-as.data.frame(tmp)
  tmp<-as.data.frame(rowMeans(tmp))
  colnames(tmp)<-i
  tmp2<-cbind(tmp2,tmp)
}  
tmp2$genes<-NULL

genes<-unique(results.tcells5$gene)
test2<-tmp2[rownames(tmp2) %in% genes,]
test2<-test2[genes,]
test2<-test2[,c(1,8,
                2,9,
                3,10,
                4,11,
                5,12,
                6,13,
                7,14)]

df<-data.frame(disease=c('Emphysema', 'Fibrosis'),
               clusters=sapply(strsplit(split='.', colnames(test2), fixed=T), '[', 2))
rownames(df)<-colnames(test2)

pheatmap(test2,
         clustering_distance_rows='correlation',
         clustering_distance_cols='correlation',
         scale='row',
         show_colnames=F,
         cluster_rows=F,
         cluster_cols=F,
         cellheight=9,
         cellwidth=16,
         gaps_col=c(2,4,6,8,10,12),
         annotation_col=df,
         breaks=seq(-2, 2, by=0.05),
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.05))))

## B cells
DefaultAssay(object=lymphoid)<-"integrated"
Idents(lymphoid)<-lymphoid$general_final_annotation
bcells<-subset(lymphoid, idents=c('B cells 1','B cells 2','B cells 3','Plasma cells','Proliferating plasma cells'))

Idents(bcells)<-'doublets'
bcells<-subset(bcells, idents='Singlet')

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
bcells<-FindNeighbors(object=bcells, reduction='pca', dims=1:15)

# Find clusters
bcells<-FindClusters(object=bcells, resolution=0.2, group.singletons=T)

# Run UMAP
bcells<-RunUMAP(object=bcells, reduction='pca', dims=1:15)
DimPlot(object=bcells, reduction='umap', label=F, pt.size=0.1)+scale_color_ucscgb()

# Find cluster markers
DefaultAssay(object=bcells)<-"RNA"
results.bcells<-all.markers(object=bcells, min.pct=0.20, log=0.25)
results.bcells<-results.bcells[results.bcells$p_val_adj<=0.05 & results.bcells$avg_log2FC>0,]
results.bcells<-results.bcells %>% arrange(cluster, -avg_log2FC)
results.bcells5<-results.bcells %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)

DotPlot(object=bcells, features=unique(results.bcells5$gene), cols='RdBu', dot.scale=6, dot.min=0.1) + RotatedAxis()

# Plot DE genes between clusters
test<-bcells@assays$RNA@data  
test<-test[rowSums(as.matrix(test))>0,]

tmp2<-data.frame(rownames(test))
colnames(tmp2)<-'genes'
bcells@meta.data$heatmap<-paste(bcells@meta.data$disease, sep='.', bcells@meta.data$integrated_snn_res.0.2)

for (i in sort(as.character(unique(bcells@meta.data$heatmap)),decreasing=F)) {
  tmp<-test[,colnames(test) %in% rownames(bcells@meta.data[bcells@meta.data$heatmap==i,])]
  tmp<-as.data.frame(tmp)
  tmp<-as.data.frame(rowMeans(tmp))
  colnames(tmp)<-i
  tmp2<-cbind(tmp2,tmp)
}  
tmp2$genes<-NULL

genes<-unique(results.bcells5$gene)
test2<-tmp2[rownames(tmp2) %in% genes,]
test2<-test2[genes,]
test2<-test2[,c(1,8,
                2,9,
                3,10,
                4,11,
                5,12,
                13,
                14,
                15,
                16)]

df<-data.frame(disease=c(rep(c('Emphysema', 'Fibrosis'), 5), rep(c('Fibrosis'), 4)),
               clusters=sapply(strsplit(split='.', colnames(test2), fixed=T), '[', 2))
rownames(df)<-colnames(test2)

pheatmap(test2,
         clustering_distance_rows='correlation',
         clustering_distance_cols='correlation',
         scale='row',
         show_colnames=F,
         cluster_rows=F,
         cluster_cols=F,
         cellheight=9,
         cellwidth=16,
         gaps_col=c(2,4,6,8,10,11,12,13),
         annotation_col=df,
         breaks=seq(-2, 2, by=0.05),
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.05))))

## NK cells
DefaultAssay(object=lymphoid)<-"integrated"
Idents(lymphoid)<-lymphoid$general_final_annotation
nks<-subset(lymphoid, idents='NK cells')

Idents(nks)<-'doublets'
nks<-subset(nks, idents='Singlet')

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
nks<-FindNeighbors(object=nks, reduction='pca', dims=1:30)

# Find clusters
nks<-FindClusters(object=nks, resolution=0.3, group.singletons=T)

# Run UMAP
nks<-RunUMAP(object=nks, reduction='pca', dims=1:30)
DimPlot(object=nks, reduction='umap', label=F, pt.size=0.1)+scale_color_aaas()

# Find cluster markers
DefaultAssay(object=nks)<-"RNA"
results.nks<-all.markers(object=nks, min.pct=0.20, log=0.25)
results.nks<-results.nks[results.nks$p_val_adj<=0.05 & results.nks$avg_log2FC>0,]
results.nks<-results.nks %>% arrange(cluster, -avg_log2FC)
results.nks5<-results.nks %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)

DotPlot(object=nks, features=unique(results.nks5$gene), cols='RdBu', dot.scale=6, dot.min=0.1) + RotatedAxis()

# Plot DE genes between clusters
test<-nks@assays$RNA@data  
test<-test[rowSums(as.matrix(test))>0,]

tmp2<-data.frame(rownames(test))
colnames(tmp2)<-'genes'
nks@meta.data$heatmap<-paste(nks@meta.data$disease, sep='.', nks@meta.data$integrated_snn_res.0.3)

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
test2<-test2[,c(1,4,
                2,5,
                3)]

df<-data.frame(disease=c(rep(c('Emphysema', 'Fibrosis'), 2),'Fibrosis'),
               clusters=c(rep('Cluster_0',each=2),rep('Cluster_1',each=2),'Cluster_2'))
rownames(df)<-colnames(test2)

pheatmap(test2,
         clustering_distance_rows='correlation',
         clustering_distance_cols='correlation',
         scale='row',
         show_colnames=F,
         cluster_rows=F,
         cluster_cols=F,
         cellheight=8,
         cellwidth=14,
         gaps_col=c(2,4),
         annotation_col=df,
         breaks=seq(-2, 2, by=0.05),
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.05))))