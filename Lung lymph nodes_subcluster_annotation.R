# Cluster 0 
DefaultAssay(object=lymphnode)<-"integrated"
Idents(lymphnode)<-lymphnode$integrated_snn_res.0.3
cluster<-subset(lymphnode, idents=0)

Idents(cluster)<-'doublets'
cluster<-subset(cluster, idents='Singlet')

# Find variable genes
cluster<-FindVariableFeatures(object=cluster, selection.method="vst", nfeatures=2000)

# Scale data
cluster<-ScaleData(object=cluster)

# Perform linear reduction    
cluster<-RunPCA(object=cluster, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=cluster, ndims=50)

# Find k-nearest neighbours
DefaultAssay(cluster)<-'integrated'
cluster<-FindNeighbors(object=cluster, reduction='pca', dims=1:30)

# Find clusters
cluster<-FindClusters(object=cluster, resolution=0.3, group.singletons=T)

# Run UMAP
cluster<-RunUMAP(object=cluster, reduction='pca', dims=1:30)
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_log2FC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis()

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='1',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='2',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='3',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='4',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='5',])),
                annotation=c(rep('T cells 1', 8627), rep('T cells_2', 1214), rep('T cells_3', 1006), rep('T cells_4', 262), rep('mixed cells', 242), rep('T cells 5', 81)),
                cluster=0)

final_annotation<-tmp

# Cluster 1
DefaultAssay(object=lymphnode)<-"integrated"
Idents(lymphnode)<-'integrated_snn_res.0.3'
cluster<-subset(lymphnode, idents=1)

Idents(cluster)<-'doublets'
cluster<-subset(cluster, idents='Singlet')

# Find variable genes
cluster<-FindVariableFeatures(object=cluster, selection.method="vst", nfeatures=2000)

# Scale data
cluster<-ScaleData(object=cluster)

# Perform linear reduction    
cluster<-RunPCA(object=cluster, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=cluster, ndims=50)

# Find k-nearest neighbours
DefaultAssay(cluster)<-'integrated'
cluster<-FindNeighbors(object=cluster, reduction='pca', dims=1:20)

# Find clusters
cluster<-FindClusters(object=cluster, resolution=0.3, group.singletons=T)

# Run UMAP
cluster<-RunUMAP(object=cluster, reduction='pca', dims=1:20)
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_log2FC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis()

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='1',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='2',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='3',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='4',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='5',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='6',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='7',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='8',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='9',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='10',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='11',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='12',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='13',])),
                annotation=c(rep('B cells_1',  4977),  rep('B cells_2', 1368), rep('B cells_3', 390), rep('B cells_4', 299), rep('B cells_5', 219), 
                             rep('B cells_6', 183), rep('B cells_7', 182), rep('B cells_8', 171), rep('B cells_9',  95), rep('B cells_10', 82), 
                             rep('B cells_11', 81), rep('B cells_12', 64), rep('B cells_13', 57), rep('B cells_14',  25)), 
                cluster=1)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 2
DefaultAssay(object=lymphnode)<-"integrated"
Idents(lymphnode)<-'integrated_snn_res.0.3'
cluster<-subset(lymphnode, idents=2)

Idents(cluster)<-'doublets'
cluster<-subset(cluster, idents='Singlet')

# Find variable genes
cluster<-FindVariableFeatures(object=cluster, selection.method="vst", nfeatures=2000)

# Scale data
cluster<-ScaleData(object=cluster)

# Perform linear reduction    
cluster<-RunPCA(object=cluster, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=cluster, ndims=50)

# Find k-nearest neighbours
DefaultAssay(cluster)<-'integrated'
cluster<-FindNeighbors(object=cluster, reduction='pca', dims=1:30)

# Find clusters
cluster<-FindClusters(object=cluster, resolution=0.3, group.singletons=T)

# Run UMAP
cluster<-RunUMAP(object=cluster, reduction='pca', dims=1:30)
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_log2FC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis()

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='1',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='2',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='3',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='4',])),
                annotation=c(rep('B cells_15', 5222), rep('B cells_16', 145), rep('mixed cells', 21), rep('B cells_17', 17), 
                             rep('B cells_18', 13)),
                cluster=2)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 3
DefaultAssay(object=lymphnode)<-"integrated"
Idents(lymphnode)<-'integrated_snn_res.0.3'
cluster<-subset(lymphnode, idents=3)

Idents(cluster)<-'doublets'
cluster<-subset(cluster, idents='Singlet')

# Find variable genes
cluster<-FindVariableFeatures(object=cluster, selection.method="vst", nfeatures=2000)

# Scale data
cluster<-ScaleData(object=cluster)

# Perform linear reduction    
cluster<-RunPCA(object=cluster, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=cluster, ndims=50)

# Find k-nearest neighbours
DefaultAssay(cluster)<-'integrated'
cluster<-FindNeighbors(object=cluster, reduction='pca', dims=1:30)

# Find clusters
cluster<-FindClusters(object=cluster, resolution=0.3, group.singletons=T)

# Run UMAP
cluster<-RunUMAP(object=cluster, reduction='pca', dims=1:30)
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_log2FC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis()

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='1',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='2',])),
                annotation=c(rep('T cells_6', 1525), rep('T cells_7', 695), rep('T cells_8', 83)),
                cluster=3)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 4
DefaultAssay(object=lymphnode)<-"integrated"
Idents(lymphnode)<-'integrated_snn_res.0.3'
cluster<-subset(lymphnode, idents=4)

Idents(cluster)<-'doublets'
cluster<-subset(cluster, idents='Singlet')

# Find variable genes
cluster<-FindVariableFeatures(object=cluster, selection.method="vst", nfeatures=2000)

# Scale data
cluster<-ScaleData(object=cluster)

# Perform linear reduction    
cluster<-RunPCA(object=cluster, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=cluster, ndims=50)

# Find k-nearest neighbours
DefaultAssay(cluster)<-'integrated'
cluster<-FindNeighbors(object=cluster, reduction='pca', dims=1:20)

# Find clusters
cluster<-FindClusters(object=cluster, resolution=0.3, group.singletons=T)

# Run UMAP
cluster<-RunUMAP(object=cluster, reduction='pca', dims=1:20)
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_log2FC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis()

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='1',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='2',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='3',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='4',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='5',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='6',])),
                annotation=c(rep('mixed cells', 894), rep('Macrophages_1', 583), rep('Macrophages_2', 291), rep('Monocytes_1', 238), 
                             rep('DCs_1', 73), rep('pDCs_1', 55), rep('Monocytes_2', 29)),
                cluster=4)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 5
DefaultAssay(object=lymphnode)<-"integrated"
Idents(lymphnode)<-'integrated_snn_res.0.3'
cluster<-subset(lymphnode, idents=5)

Idents(cluster)<-'doublets'
cluster<-subset(cluster, idents='Singlet')

# Find variable genes
cluster<-FindVariableFeatures(object=cluster, selection.method="vst", nfeatures=2000)

# Scale data
cluster<-ScaleData(object=cluster)

# Perform linear reduction    
cluster<-RunPCA(object=cluster, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=cluster, ndims=50)

# Find k-nearest neighbours
DefaultAssay(cluster)<-'integrated'
cluster<-FindNeighbors(object=cluster, reduction='pca', dims=1:30)

# Find clusters
cluster<-FindClusters(object=cluster, resolution=0.3, group.singletons=T)

# Run UMAP
cluster<-RunUMAP(object=cluster, reduction='pca', dims=1:30)
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_log2FC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis()

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='1',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='2',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='3',])),
                annotation=c(rep('NK cells_1', 1618), rep('NK cells_2', 355), rep('mixed cells', 151), rep('NK cells_3', 11)),
                cluster=5)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 6
DefaultAssay(object=lymphnode)<-"integrated"
Idents(lymphnode)<-'integrated_snn_res.0.3'
cluster<-subset(lymphnode, idents=6)

Idents(cluster)<-'doublets'
cluster<-subset(cluster, idents='Singlet')

# Find variable genes
cluster<-FindVariableFeatures(object=cluster, selection.method="vst", nfeatures=2000)

# Scale data
cluster<-ScaleData(object=cluster)

# Perform linear reduction    
cluster<-RunPCA(object=cluster, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=cluster, ndims=50)

# Find k-nearest neighbours
DefaultAssay(cluster)<-'integrated'
cluster<-FindNeighbors(object=cluster, reduction='pca', dims=1:20)

# Find clusters
cluster<-FindClusters(object=cluster, resolution=0.3, group.singletons=T)

# Run UMAP
cluster<-RunUMAP(object=cluster, reduction='pca', dims=1:20)
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_log2FC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis()

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='1',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='2',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='3',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='4',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='5',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='6',])),
                annotation=c(rep('Plasma cells_1', 706), rep('Plasma cells_2', 289), rep('Plasma cells_3', 196), rep('Plasma cells_4', 169), 
                             rep('Plasma cells_5', 109), rep('Plasma cells_6', 104), rep('Plasma cells_7', 90)), 
                cluster=6)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 7
DefaultAssay(object=lymphnode)<-"integrated"
Idents(lymphnode)<-'integrated_snn_res.0.3'
cluster<-subset(lymphnode, idents=7)

Idents(cluster)<-'doublets'
cluster<-subset(cluster, idents='Singlet')

# Find variable genes
cluster<-FindVariableFeatures(object=cluster, selection.method="vst", nfeatures=2000)

# Scale data
cluster<-ScaleData(object=cluster)

# Perform linear reduction    
cluster<-RunPCA(object=cluster, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=cluster, ndims=50)

# Find k-nearest neighbours
DefaultAssay(cluster)<-'integrated'
cluster<-FindNeighbors(object=cluster, reduction='pca', dims=1:20)

# Find clusters
cluster<-FindClusters(object=cluster, resolution=0.3, group.singletons=T)

# Run UMAP
cluster<-RunUMAP(object=cluster, reduction='pca', dims=1:20)
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_log2FC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis()

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='1',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='2',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='3',])),
                annotation=c(rep('Proliferating T cells_1', 258), rep('Proliferating T/B cells_1', 199), rep('Proliferating plasma cells_1', 98), 
                             rep('Proliferating macrophages_1', 34)), 
                cluster=7)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 8
DefaultAssay(object=lymphnode)<-"integrated"
Idents(lymphnode)<-'integrated_snn_res.0.3'
cluster<-subset(lymphnode,  idents=8)

Idents(cluster)<-'doublets'
cluster<-subset(cluster, idents='Singlet')

# Find variable genes
cluster<-FindVariableFeatures(object=cluster, selection.method="vst", nfeatures=2000)

# Scale data
cluster<-ScaleData(object=cluster)

# Perform linear reduction    
cluster<-RunPCA(object=cluster, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=cluster, ndims=50)

# Find k-nearest neighbours
DefaultAssay(cluster)<-'integrated'
cluster<-FindNeighbors(object=cluster, reduction='pca', dims=1:30)

# Find clusters
cluster<-FindClusters(object=cluster, resolution=0.3, group.singletons=T)

# Run UMAP
cluster<-RunUMAP(object=cluster, reduction='pca', dims=1:30)
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_log2FC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis()

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='1',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='2',])),
                annotation=c(rep('Mast cells_1', 301), rep('Mast cells_2', 190), rep('mixed cells', 19)),
                cluster=8)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 9
DefaultAssay(object=lymphnode)<-"integrated"
Idents(lymphnode)<-'integrated_snn_res.0.3'
cluster<-subset(lymphnode, idents=9)

Idents(cluster)<-'doublets'
cluster<-subset(cluster, idents='Singlet')

# Find variable genes
cluster<-FindVariableFeatures(object=cluster, selection.method="vst", nfeatures=2000)

# Scale data
cluster<-ScaleData(object=cluster)

# Perform linear reduction    
cluster<-RunPCA(object=cluster, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=cluster, ndims=50)

# Find k-nearest neighbours
DefaultAssay(cluster)<-'integrated'
cluster<-FindNeighbors(object=cluster, reduction='pca', dims=1:30)

# Find clusters
cluster<-FindClusters(object=cluster, resolution=0.3, group.singletons=T)

# Run UMAP
cluster<-RunUMAP(object=cluster, reduction='pca', dims=1:30)
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_log2FC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis()

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',])),
                annotation=c(rep('B cells_19', 282)),
                cluster=9)

final_annotation<-rbind(final_annotation, tmp)

# Save the annotation in the dataset
test<-final_annotation
test$annotation<-as.character(test$annotation)
test$merge<-paste0(test$annotation, '_cluster ', test$cluster)

lymphnode$detailed_final_annotation<-''
lymphnode$detailed_final_annotation<-test$merge[match(rownames(lymphnode@meta.data), test$barcodes)]
lymphnode@meta.data[lymphnode@meta.data$doublets=='Doublet',]$detailed_final_annotation<-'Doublet'

# Reduce to major annotations
lymphnode$general_final_annotation<-''
lymphnode@meta.data[lymphnode@meta.data$doublets=='Doublet',]$general_final_annotation<-'Doublet'

table(lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==0,]$detailed_final_annotation)
lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==0 & grepl('T cells', lymphnode@meta.data$detailed_final_annotation),]$general_final_annotation<-'T cells 1'
lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==0 & grepl('mixed', lymphnode@meta.data$detailed_final_annotation),]$general_final_annotation<-'mixed cells'

table(lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==1,]$detailed_final_annotation)
lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==1 & grepl('B cells', lymphnode@meta.data$detailed_final_annotation),]$general_final_annotation<-'B cells 1'

table(lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==2,]$detailed_final_annotation)
lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==2 & grepl('B cells', lymphnode@meta.data$detailed_final_annotation),]$general_final_annotation<-'B cells 2'
lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==2 & grepl('mixed', lymphnode@meta.data$detailed_final_annotation),]$general_final_annotation<-'mixed cells'

table(lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==3,]$detailed_final_annotation)
lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==3 & grepl('T cells', lymphnode@meta.data$detailed_final_annotation),]$general_final_annotation<-'T cells 2'

table(lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==4,]$detailed_final_annotation)
lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==4 & grepl('Macrophages', lymphnode@meta.data$detailed_final_annotation),]$general_final_annotation<-'Macrophages'
lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==4 & grepl('Monocytes', lymphnode@meta.data$detailed_final_annotation),]$general_final_annotation<-'Monocytes'
lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==4 & grepl('DCs', lymphnode@meta.data$detailed_final_annotation),]$general_final_annotation<-'DCs'
lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==4 & grepl('pDCs', lymphnode@meta.data$detailed_final_annotation),]$general_final_annotation<-'plasmacytoid DCs'
lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==4 & grepl('mixed', lymphnode@meta.data$detailed_final_annotation),]$general_final_annotation<-'mixed cells'

table(lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==5,]$detailed_final_annotation)
lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==5 & grepl('NK cells', lymphnode@meta.data$detailed_final_annotation),]$general_final_annotation<-'NK cells'
lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==5 & grepl('mixed', lymphnode@meta.data$detailed_final_annotation),]$general_final_annotation<-'mixed cells'

table(lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==6,]$detailed_final_annotation)
lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==6 & grepl('Plasma cells', lymphnode@meta.data$detailed_final_annotation),]$general_final_annotation<-'Plasma cells'

table(lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==7,]$detailed_final_annotation)
lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==7 & grepl('T/B cells', lymphnode@meta.data$detailed_final_annotation),]$general_final_annotation<-'Proliferating T/B cells'
lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==7 & grepl('macrophages', lymphnode@meta.data$detailed_final_annotation),]$general_final_annotation<-'Proliferating macrophages'
lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==7 & grepl('plasma cells', lymphnode@meta.data$detailed_final_annotation),]$general_final_annotation<-'Proliferating plasma cells'
lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==7 & grepl('T cells', lymphnode@meta.data$detailed_final_annotation),]$general_final_annotation<-'Proliferating T cells'

table(lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==8,]$detailed_final_annotation)
lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==8 & grepl('Mast cells', lymphnode@meta.data$detailed_final_annotation),]$general_final_annotation<-'Mast cells'
lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==8 & grepl('mixed', lymphnode@meta.data$detailed_final_annotation),]$general_final_annotation<-'mixed cells'

table(lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==9,]$detailed_final_annotation)
lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==9 & lymphnode@meta.data$detailed_final_annotation!='Doublet',]$general_final_annotation<-'B cells 3'

# Check final results and refine
for (i in c(0:9)) {
  print(paste0('Cluster ', i))
  print(table(lymphnode@meta.data[lymphnode@meta.data$integrated_snn_res.0.3==i,]$general_final_annotation))
}

lymphnode$general_final_annotation<-factor(lymphnode$general_final_annotation)
lymphnode$general_final_annotation<-factor(lymphnode$general_final_annotation, 
                                               levels=levels(lymphnode$general_final_annotation)[c(6,9,4,12,
                                                                                                   7,
                                                                                                   17,18,
                                                                                                   1,2,3,11,
                                                                                                   10,
                                                                                                   13,15,16,14,
                                                                                                   5,8)])


DimPlot(lymphnode, group.by='general_final_annotation') + scale_color_manual(values=c("#8B4513","#D2691E","#CD853F","#F4A460","#CD5C5C","#228B22","#32CD32",
                                                                                      "#000080","#0000FF","#4682B4","#1E90FF","#2F4F4F","#9370DB",
                                                                                      "#9932CC","#FF00FF","#D02090","white","white"))

# Calculate patient proportions
set.seed(12345)
downsampling<-lymphnode@meta.data %>% group_by(LuT) %>% sample_n(305, replace=F)
data<-round(prop.table(x=table(downsampling$LuT, downsampling$general_final_annotation), margin=1), digits=3)
data<-melt(data)
data$Condition<-lymhnode$disease[match(data$Var1, lymphnode$LuT)]
data<-data[data$Condition %in% c('Fibrosis', 'Emphysema'),]
data$Condition<-factor(data$Condition, levels=c('Emphysema', 'Fibrosis'))
data<-data[data$Var2 %in% levels(data$Var2)[1:17],]

ggplot(data, aes(x=Condition, y=value, fill=Condition))+
  stat_boxplot(geom='errorbar', linetype=1, width=0.8)+ 
  geom_boxplot(outlier.shape=1)+   
  ylab('')+
  xlab('')+
  facet_wrap(~Var2, ncol=6, scale='free')+
  scale_fill_manual(values=c("#1B9E77", "#7570B3"))+
  theme(axis.text.x=element_text(size=14, color='black'),
        axis.text.y=element_text(size=14, color='black'),
        strip.text=element_text(size=12),
        strip.background=element_rect(fill='gray90', colour='black'),
        legend.position='top',
        legend.key=element_rect(colour='black', fill='white'),
        legend.spacing.x=unit(0.3, 'cm'),
        panel.background=element_rect(fill='white'),
        panel.border=element_blank(),
        axis.line=element_line(colour="black"))

# Perform statistical analysis 
# Perform statistical analysis for cluster proportions
for (i in levels(data$Var2)[c(1,2,5:8,10,11,13:16)]){
  test<-data[data$Var2==i,]
  print(paste0('Significance test for subset ',i))
  a<-shapiro.test(test[test$Condition=='Emphysema',]$value)
  b<-shapiro.test(test[test$Condition=='Fibrosis',]$value)
  if(a$p.value <0.05 & b$p.value >=0.05){
    print('Emphysema values failed normality test')
    d<-wilcox.test(x=test[test$Condition=='Emphysema',]$value, y=test[test$Condition=='Fibrosis',]$value)  
  }else if(a$p.value >=0.05 & b$p.value <0.05){
    print('Fibrosis values failed normality test')
    d<-wilcox.test(x=test[test$Condition=='Emphysema',]$value, y=test[test$Condition=='Fibrosis',]$value) 
  }else if(a$p.value <0.05 & b$p.value <0.05){
    print('Both Emphysema and Fibrosis values failed normality test')
    d<-wilcox.test(x=test[test$Condition=='Emphysema',]$value, y=test[test$Condition=='Fibrosis',]$value)
  }else{
    print('Both Emphysema and Fibrosis values come from a normal distribution')
    d<-t.test(x=test[test$Condition=='Emphysema',]$value, y=test[test$Condition=='Fibrosis',]$value)
  }
  print(d)
}

# Calculate DE gene markers
DefaultAssay(object=lymphnode)<-"RNA"
Idents(lymphnode)<-'general_final_annotation'

results<-all.markers(object=lymphnode, min.pct=0.20, log=0.4)
results<-results[results$p_val_adj <= 0.05 & results$avg_log2FC>0,]
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)
results5$cluster<-factor(results5$cluster, levels=levels(results5$cluster)[c(6,9,4,12,7,17,18,1,2,3,11,10,13,
                                                                             15,16,14,5,8)])

DotPlot(object=lymphnode, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis()

test<-lymphnode@assays$RNA@data  
tmp2<-data.frame(rownames(test))
colnames(tmp2)<-'genes'

lymphnode$heatmap<-paste(lymphnode$general_final_annotation, sep='_', lymphnode$disease)
for (i in unique(lymphnode@meta.data$heatmap)) {
  tmp<-test[,colnames(test) %in% rownames(lymphnode@meta.data[lymphnode@meta.data$heatmap==i,])]
  tmp<-as.data.frame(tmp)
  tmp<-as.data.frame(rowMeans(tmp))
  colnames(tmp)<-i
  tmp2<-cbind(tmp2,tmp)
}  
tmp2$genes<-NULL

test2<-tmp2[rownames(tmp2) %in% unique(results5$gene[1:80]),]
test2<-test2[unique(results5$gene[1:80]),]
test2<-test2[,c(11,29,
                14,30,
                15,33,
                13,32,
                5,28,
                1,17,
                7,23,
                4,22,
                3,21,
                35,36,
                8,25,
                2,20,
                34,31,
                9,24,
                12,19,
                16,26)]

df<-data.frame(celltype=sapply(strsplit(colnames(test2), '_'), '[', 1))   
rownames(df)<-colnames(test2)

anno<-c("#8B4513","#D2691E","#CD853F","#F4A460","#CD5C5C","#228B22","#32CD32","#000080","#0000FF","#4682B4",
        "#1E90FF","#2F4F4F","#9370DB","#9932CC","#FF00FF","#D02090")
names(anno)<-unique(df$celltype)

pheatmap(test2,
         clustering_distance_rows='correlation',
         clustering_distance_cols='correlation',
         scale='row',
         show_rownames=T,
         show_colnames=F,
         cluster_rows=F,
         cluster_cols=F,
         breaks=seq(-2, 2, by=0.1),
         gaps_col=c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30),
         cellheight=6,
         cellwidth=9,
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.1))),
         annotation_col=df,
         annotation_colors=list(celltype=anno))

# Add clinical information
df<-read.delim('261020_Hannover_sample_table.txt', header=T)
df<-df[,c(1,2,4:11,22,24,27:32)]

lymphnode[['disease']]<-df$disease[match(lymphnode@meta.data$LuT, df$LuT)]
lymphnode[['disease_tissue']]<-df$disease_tissue[match(lymphnode@meta.data$LuT, df$LuT)]
lymphnode[['classification']]<-df$Classification[match(lymphnode@meta.data$LuT, df$LuT)]
lymphnode[['other_diagnosis']]<-df$other_diagnosis[match(lymphnode@meta.data$LuT, df$LuT)]
lymphnode[['site']]<-df$site[match(lymphnode@meta.data$LuT, df$LuT)]
lymphnode[['immediate']]<-df$Immediate[match(lymphnode@meta.data$LuT, df$LuT)]
lymphnode[['age']]<-df$Age[match(lymphnode@meta.data$LuT, df$LuT)]
lymphnode[['sex']]<-df$sex[match(lymphnode@meta.data$LuT, df$LuT)]
lymphnode[['BMI']]<-df$BMI[match(lymphnode@meta.data$LuT, df$LuT)]
lymphnode[['blood_group']]<-df$Blood_group[match(lymphnode@meta.data$LuT, df$LuT)]
lymphnode[['CMV_vaccination']]<-df$IgG_CMV_Empf[match(lymphnode@meta.data$LuT, df$LuT)]
lymphnode[['FEV1']]<-df$FEV1[match(lymphnode@meta.data$LuT, df$LuT)]
lymphnode[['FVC']]<-df$FVC1[match(lymphnode@meta.data$LuT, df$LuT)]
lymphnode[['smoking']]<-df$Nicotine_abuse[match(lymphnode@meta.data$LuT, df$LuT)]