# Cluster 0 
DefaultAssay(object=combined_list)<-"integrated"
Idents(combined_list)<-combined_list$integrated_snn_res.0.4
cluster<-subset(combined_list, idents=0)

Idents(cluster)<-'doublets_individual'
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
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
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
                annotation=c(rep('T cells 1', 4680), rep('T cells_2', 3825), rep('NK cells', 3095), rep('T cells_3', 1076),
                             rep('T cells_4', 201), rep('T cells_5', 97), rep('T cells_6', 20)),
                cluster=0)

final_annotation<-tmp

# Cluster 1
DefaultAssay(object=combined_list)<-"integrated"
Idents(combined_list)<-'integrated_snn_res.0.4'
cluster<-subset(combined_list, idents=1)

Idents(cluster)<-'doublets_individual'
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
cluster<-FindNeighbors(object=cluster, reduction='pca', dims=1:15)

# Find clusters
cluster<-FindClusters(object=cluster, resolution=0.3, group.singletons=T)

# Run UMAP
cluster<-RunUMAP(object=cluster, reduction='pca', dims=1:15)
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
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
                annotation=c(rep('Macrophages_1', 4060), rep('Macrophages_2', 3205), rep('Macrophages_3', 1401),
                             rep('Macrophages_4', 502), rep('Macrophages_5', 19)),
                cluster=1)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 2 
DefaultAssay(object=combined_list)<-"integrated"
Idents(combined_list)<-'integrated_snn_res.0.4'
cluster<-subset(combined_list, idents=2)

Idents(cluster)<-'doublets_individual'
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
cluster<-FindNeighbors(object=cluster, reduction='pca', dims=1:10)

# Find clusters
cluster<-FindClusters(object=cluster, resolution=0.3, group.singletons=T)

# Run UMAP
cluster<-RunUMAP(object=cluster, reduction='pca', dims=1:10)
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8 ,dot.min=0.1) + RotatedAxis()

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='1',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='2',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='3',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='4',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='5',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='6',])),
                annotation=c(rep('mixed', 3381), rep('Plasma cells_1',  1999), rep('Plasma cells_2', 411), rep('Plasma cells_3', 219),
                             rep('Plasma cells_4', 82), rep('Plasma cells_5', 44), rep('mixed', 10)),
                cluster=2)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 3
DefaultAssay(object=combined_list)<-"integrated"
Idents(combined_list)<-'integrated_snn_res.0.4'
cluster<-subset(combined_list, idents=3)

Idents(cluster)<-'doublets_individual'
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
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis()

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='1',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='2',])),
                annotation=c(rep('Macrophages_6', 1764), rep('Macrophages_7', 1061), rep('Macrophages_8', 613)),
                cluster=3)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 4
DefaultAssay(object=combined_list)<-"integrated"
Idents(combined_list)<-'integrated_snn_res.0.4'
cluster<-subset(combined_list, idents=4)

Idents(cluster)<-'doublets_individual'
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
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
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
                annotation=c(rep('Monocytes_1', 1962), rep('Monocytes_2', 1094), rep('Monocytes_3', 938),
                             rep('Monocytes_4', 763), rep('Monocytes_5', 51)),
                cluster=4)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 5 
DefaultAssay(object=combined_list)<-"integrated"
Idents(combined_list)<-'integrated_snn_res.0.4'
cluster<-subset(combined_list, idents=5)

Idents(cluster)<-'doublets_individual'
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
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
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
                annotation=c(rep('mixed', 1045), rep('Club cells_1', 869), rep('Club cells_2', 655),
                             rep('Club cells 3', 520), rep('AT1 cells', 189), rep('Club cells_4', 46)),
                cluster=5)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 6 
DefaultAssay(object=combined_list)<-"integrated"
Idents(combined_list)<-'integrated_snn_res.0.4'
cluster<-subset(combined_list, idents=6)

Idents(cluster)<-'doublets_individual'
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
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8,dot.min=0.1) + RotatedAxis()

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
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='8',])),
                annotation=c(rep('Ciliated cells_1', 640), rep('mixed', 626), rep('mixed', 475), rep('Ciliated cells_2', 426),
                             rep('Ciliated cells_3', 342), rep('Ciliated cells_4', 153), rep('mixed', 147), rep('Ciliated cells_5', 37)),
                cluster=6)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 7
DefaultAssay(object=combined_list)<-"integrated"
Idents(combined_list)<-'integrated_snn_res.0.4'
cluster<-subset(combined_list, idents=7)

Idents(cluster)<-'doublets_individual'
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
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
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
                annotation=c(rep('Neutrophils_1', 1273), rep('Neutrophils_2', 772), rep('Neutrophils_3', 442), rep('Eosinophils', 105),
                             rep('Neutrophils_4', 30)),
                cluster=7)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 8
DefaultAssay(object=combined_list)<-"integrated"
Idents(combined_list)<-'integrated_snn_res.0.4'
cluster<-subset(combined_list, idents=8)

Idents(cluster)<-'doublets_individual'
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
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8 ,dot.min=0.1) + RotatedAxis()

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='1',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='2',])),
                annotation=c(rep('Macrophages_9', 1135), rep('Macrophages_10', 999), rep('Macrophages_11', 49)),
                cluster=8)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 9
DefaultAssay(object=combined_list)<-"integrated"
Idents(combined_list)<-'integrated_snn_res.0.4'
cluster<-subset(combined_list, idents=9)

Idents(cluster)<-'doublets_individual'
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
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis()

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='1',])),
                annotation=c(rep('mixed', 1162), rep('T cells_7', 867)),
                cluster=9)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 10
DefaultAssay(object=combined_list)<-"integrated"
Idents(combined_list)<-'integrated_snn_res.0.4'
cluster<-subset(combined_list, idents=10)

Idents(cluster)<-'doublets_individual'
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
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8,dot.min=0.1) + RotatedAxis()

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='1',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='2',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='3',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='4',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='5',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='6',])),
                annotation=c(rep('Erythrocytes_1', 1280), rep('Erythrocytes_2', 437), rep('Erythrocytes_3', 99),
                             rep('Erythrocytes_4', 40), rep('mixed', 32), rep('Erythrocytes_5', 26),
                             rep('Erythrocytes_6', 21)),
                cluster=10)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 11
DefaultAssay(object=combined_list)<-"integrated"
Idents(combined_list)<-'integrated_snn_res.0.4'
cluster<-subset(combined_list, idents=11)

Idents(cluster)<-'doublets_individual'
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
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis()

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='1',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='2',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='3',])),
                annotation=c(rep('Dendritic cells_1', 666), rep('Dendritic cells_2', 627), rep('Dendritic cells_3', 405),
                             rep('Dendritic cells_4', 42)),
                cluster=11)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 12
DefaultAssay(object=combined_list)<-"integrated"
Idents(combined_list)<-'integrated_snn_res.0.4'
cluster<-subset(combined_list, idents=12)

Idents(cluster)<-'doublets_individual'
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
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='1',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='2',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='3',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='4',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='5',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='6',])),
                annotation=c(rep('AT2 cells_1', 619), rep('AT2 cells_2', 464), rep('AT2 cells_3', 95) ,rep('AT2 cells_4', 84),
                             rep('AT2 cells_5', 59), rep('mixed', 24), rep('mixed', 19)),
                cluster=12)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 13
DefaultAssay(object=combined_list)<-"integrated"
Idents(combined_list)<-'integrated_snn_res.0.4'
cluster<-subset(combined_list, idents=13)

Idents(cluster)<-'doublets_individual'
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
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
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
                annotation=c(rep('Cycling T/NK cells', 370), rep('Cycling macrophages_1', 339), rep('Cycling macrophages_2', 280),
                             rep('Cycling dendritic cells', 28), rep('Cycling macrophages_3', 18)),
                cluster=13)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 14
DefaultAssay(object=combined_list)<-"integrated"
Idents(combined_list)<-'integrated_snn_res.0.4'
cluster<-subset(combined_list, idents=14)

Idents(cluster)<-'doublets_individual'
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
cluster<-FindNeighbors(object=cluster, reduction='pca', dims=1:40)

# Find clusters
cluster<-FindClusters(object=cluster, resolution=0.3, group.singletons=T)

# Run UMAP
cluster<-RunUMAP(object=cluster, reduction='pca', dims=1:40)
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis()

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='1',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='2',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='3',])),
                annotation=c(rep('B cells_1', 559), rep('B cells_2', 362), rep('B cells_3', 321), rep('B cells_4', 15)),
                cluster=14)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 15
DefaultAssay(object=combined_list)<-"integrated"
Idents(combined_list)<-'integrated_snn_res.0.4'
cluster<-subset(combined_list, idents=15)

Idents(cluster)<-'doublets_individual'
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
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8,dot.min=0.1) + RotatedAxis()

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='1',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='2',])),
                annotation=c(rep('CD56 NK cells_1', 746), rep('CD56 NK cells_2', 402), rep('CD56 NK cells_3', 27)),
                cluster=15)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 16
DefaultAssay(object=combined_list)<-"integrated"
Idents(combined_list)<-'integrated_snn_res.0.4'
cluster<-subset(combined_list, idents=16)

Idents(cluster)<-'doublets_individual'
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
cluster<-FindNeighbors(object=cluster, reduction='pca', dims=1:40)

# Find clusters
cluster<-FindClusters(object=cluster, resolution=0.3, group.singletons=T)

# Run UMAP
cluster<-RunUMAP(object=cluster, reduction='pca', dims=1:40)
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis()

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='1',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='2',])),
                annotation=c(rep('Mast cells/Basophils_1', 427), rep('Mast cells/Basophils_2', 199), rep('Mast cells/Basophils_3', 182)),
                cluster=16)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 17
DefaultAssay(object=combined_list)<-"integrated"
Idents(combined_list)<-'integrated_snn_res.0.4'
cluster<-subset(combined_list, idents=17)

Idents(cluster)<-'doublets_individual'
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
cluster<-FindNeighbors(object=cluster, reduction='pca', dims=1:40)

# Find clusters
cluster<-FindClusters(object=cluster, resolution=0.3, group.singletons=T)

# Run UMAP
cluster<-RunUMAP(object=cluster, reduction='pca', dims=1:40)
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8,dot.min=0.1) + RotatedAxis()

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='1',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='2',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='3',])),
                annotation=c(rep('Immature neutrophils_1', 153), rep('Immature neutrophils_2', 97), rep('Immature neutrophils_3', 75),
                             rep('Immature neutrophils_4', 56)),
                cluster=17)

final_annotation<-rbind(final_annotation, tmp)

# Cluster 18
DefaultAssay(object=combined_list)<-"integrated"
Idents(combined_list)<-'integrated_snn_res.0.4'
cluster<-subset(combined_list, idents=18)

Idents(cluster)<-'doublets_individual'
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
cluster<-FindNeighbors(object=cluster, reduction='pca', dims=1:40)

# Find clusters
cluster<-FindClusters(object=cluster, resolution=0.3, group.singletons=T)

# Run UMAP
cluster<-RunUMAP(object=cluster, reduction='pca', dims=1:40)
DimPlot(object=cluster, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=cluster)<-"RNA"
results<-all.markers(object=cluster, min.pct=0.25, log=0.4)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
results<-results %>% arrange(cluster, -pct.1)
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=pct.1)

DotPlot(object=cluster, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis()

# Save barcode annotation
table(cluster$integrated_snn_res.0.3)

tmp<-data.frame(barcodes=c(rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='0',]),
                           rownames(cluster@meta.data[cluster@meta.data$integrated_snn_res.0.3=='1',])),
                annotation=c(rep('Lipofibroblasts', 110), rep('Matrix fibroblasts', 57)),
                cluster=18)

final_annotation<-rbind(final_annotation, tmp)

# Save the annotation in the dataset
test<-final_annotation
test$annotation<-as.character(test$annotation)
test$merge<-paste0(test$annotation, '_cluster ', test$cluster)

combined_list$detailed_final_annotation<-''
combined_list$detailed_final_annotation<-test$merge[match(rownames(combined_list@meta.data),test$barcodes)]
combined_list@meta.data[combined_list@meta.data$doublets_individual=='Doublet',]$detailed_final_annotation<-'Doublet'

# Reduce to major annotations
combined_list$general_final_annotation<-''
combined_list@meta.data[combined_list@meta.data$doublets_individual=='Doublet',]$general_final_annotation<-'Doublet'

table(combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==0,]$detailed_final_annotation)
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==0 & grepl('T cells', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'T cells 1'
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==0 & grepl('NK cells', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'NK cells'

table(combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==1,]$detailed_final_annotation)
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==1 & combined_list@meta.data$detailed_final_annotation!='Doublet',]$general_final_annotation<-'Macrophages 1'

table(combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==2,]$detailed_final_annotation)
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==2 & grepl('Plasma', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'Plasma cells'
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==2 & grepl('mixed', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'mixed cells'

table(combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==3,]$detailed_final_annotation)
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==3 & combined_list@meta.data$detailed_final_annotation!='Doublet',]$general_final_annotation<-'Macrophages 2'

table(combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==4,]$detailed_final_annotation)
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==4 & combined_list@meta.data$detailed_final_annotation!='Doublet',]$general_final_annotation<-'Monocytes'

table(combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==5,]$detailed_final_annotation)
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==5 & grepl('AT1', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'AT1 cells'
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==5 & grepl('Club', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'Club cells'
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==5 & grepl('mixed', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'mixed cells'

table(combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==6,]$detailed_final_annotation)
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==6 & grepl('Ciliated', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'Ciliated cells'
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==6 & grepl('mixed', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'mixed cells'

table(combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==7,]$detailed_final_annotation)
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==7 & grepl('Neutrophils', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'Neutrophils'
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==7 & grepl('Eosinophils', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'Eosinophils'

table(combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==8,]$detailed_final_annotation)
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==8 & combined_list@meta.data$detailed_final_annotation!='Doublet',]$general_final_annotation<-'Macrophages 3'

table(combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==9,]$detailed_final_annotation)
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==9 & grepl('T cells', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'T cells 2'
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==9 & grepl('mixed', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'mixed cells'

table(combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==10,]$detailed_final_annotation)
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==10 & grepl('Erythrocytes', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'Erythrocytes'
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==10 & grepl('mixed', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'mixed cells'

table(combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==11,]$detailed_final_annotation)
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==11 & grepl('Dendritic', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'Dendritic cells'

table(combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==12,]$detailed_final_annotation)
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==12 & grepl('AT2', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'AT2 cells'
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==12 & grepl('mixed', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'mixed cells'

table(combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==13,]$detailed_final_annotation)
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==13 & grepl('dendritic', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'Cycling dendritic cells'
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==13 & grepl('macrophages', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'Cycling macrophages'
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==13 & grepl('T/NK', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'Cycling T/NK cells'

table(combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==14,]$detailed_final_annotation)
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==14 & combined_list@meta.data$detailed_final_annotation!='Doublet',]$general_final_annotation<-'B cells'

table(combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==15,]$detailed_final_annotation)
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==15 & combined_list@meta.data$detailed_final_annotation!='Doublet',]$general_final_annotation<-'CD56 NK cells'

table(combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==16,]$detailed_final_annotation)
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==16 & combined_list@meta.data$detailed_final_annotation!='Doublet',]$general_final_annotation<-'Mast cells/Basophils'

table(combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==17,]$detailed_final_annotation)
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==17 & combined_list@meta.data$detailed_final_annotation!='Doublet',]$general_final_annotation<-'Immature neutrophils'

table(combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==18,]$detailed_final_annotation)
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==18 & grepl('Alveolar', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'Alveolar fibroblasts'
combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==18 & grepl('Adventitial', combined_list@meta.data$detailed_final_annotation),]$general_final_annotation<-'Adventitial fibroblasts'

# Check final results and refine
for (i in c(0:18)) {
  print(paste0('Cluster ', i))
  print(table(combined_list@meta.data[combined_list@meta.data$integrated_snn_res.0.4==i,]$general_final_annotation))
}

combined_list$general_final_annotation<-factor(combined_list$general_final_annotation)
combined_list$general_final_annotation<-factor(combined_list$general_final_annotation, 
                                               levels=levels(combined_list$general_final_annotation)[c(16,17,18,22,10,23,14,12,19,26,27,3,25,24,4,8,7,9,1,2,6,5,13,15,20,11,21)])

DimPlot(combined_list, group.by='general_final_annotation') + scale_color_manual(values=c("#DBDB29","#BCBD22","#82853B","#FFC156","#FFB579","#FE870D","#FDA440",
                                                                                          "#ED8F47","#B15928","#ED4F50","#E4201F","#B294C7","#865FAB","#D4A55B",
                                                                                          "#B89B74","#99CD91","#69B764","#309343","#A6CEE3","#68A6CD","#569EA4",
                                                                                          "#2A7FB7","#595959","#BD0A36","#FFB6B0","white","white"))

# Calculate patient proportions
set.seed(12345)
downsampling<-combined_list@meta.data %>% group_by(LuT) %>% sample_n(243, replace=F)
data<-round(prop.table(x=table(downsampling$LuT, downsampling$general_final_annotation), margin=1), digits=3)
data<-melt(data)
data$Condition<-combined_list$disease[match(data$Var1, combined_list$LuT)]
data<-data[data$Condition %in% c('Fibrosis', 'Emphysema', 'Tu.free'),]
data$Condition<-factor(data$Condition, levels=c('Tu.free', 'Emphysema', 'Fibrosis'))
data<-data[data$Var2 %in% levels(data$Var2)[19:25],]

ggplot(data, aes(x=Condition,y=value,fill=Condition))+
  stat_boxplot(geom='errorbar', linetype=1, width=0.8)+ 
  geom_boxplot(outlier.shape=1)+   
  ylab('')+
  xlab('')+
  facet_wrap(~Var2, ncol=5, scale='free')+
  scale_fill_manual(values=c("#E6AB02", "#1B9E77", "#7570B3"))+
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
for (i in levels(data$Var2)){
  df<-data[data$Var2==i,]
  print(paste0('Test result for ',i))
  res.aov<-aov(value~Condition, df)
  var<-leveneTest(value~Condition, df)
  var<-var$`Pr(>F)`[1]
  resi<-shapiro.test(residuals(object=res.aov))
  resi<-resi$p.value
  if(var | resi <= 0.05){
    dunn.test(x=df$value,g=df$Condition)
  }else{
    print(TukeyHSD(res.aov))
  }
}

# Calculate DE gene markers
DefaultAssay(object=combined_list)<-"RNA"

Idents(combined_list)<-'general_final_annotation'

results<-all.markers(object=combined_list, min.pct=0.20, log=0.4)
results<-results[results$p_val_adj <= 0.05 & results$avg_log2FC>0,]
results3<-results %>% group_by(cluster) %>% top_n(n=3, wt=avg_log2FC)
results3$cluster<-factor(results3$cluster, levels=levels(results3$cluster)[c(1,3,2,8,7,12,27,26,9,19,23,21,16,18,15,6,24,11,25,5,13,10,17,20,22,14,4)])

DotPlot(combined_list, features=unique(results5$gene), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis()

test<-combined_list@assays$RNA@data  
tmp2<-data.frame(rownames(test))
colnames(tmp2)<-'genes'

combined_list$heatmap<-paste(combined_list$general_final_annotation, sep='_', combined_list$disease)
for (i in unique(combined_list@meta.data$heatmap)) {
  tmp<-test[,colnames(test) %in% rownames(combined_list@meta.data[combined_list@meta.data$heatmap==i,])]
  tmp<-as.data.frame(tmp)
  tmp<-as.data.frame(rowMeans(tmp))
  colnames(tmp)<-i
  tmp2<-cbind(tmp2,tmp)
}  
tmp2$genes<-NULL

test2<-tmp2[rownames(tmp2) %in% unique(results3 %>% arrange(cluster) %>% pull(gene)),]
test2<-test2[unique(results3 %>% arrange(cluster) %>% pull(gene)),]
test2<-test2[,c(94,1,64,27,
                87,3,75,42,
                85,2,69,33,
                98,8,57,29,
                99,7,72,28,
                96,12,70,31,
                104,54,79,48,
                107,26,67,44,
                100,9,63,32,
                86,19,55,38,
                91,23,74,51,
                102,21,73,41,
                89,16,61,34,
                93,18,65,40,
                95,15,59,47,
                90,6,68,39,
                106,24,81,52,
                101,11,66,45,
                103,25,77,49,
                97,5,76,30,             
                84,13,71,35,
                83,10,60,53,
                92,17,58,43,
                105,20,78,50,
                22,80,46)]


df<-data.frame(celltype=sapply(strsplit(colnames(test2), '_'), '[', 1))   
rownames(df)<-colnames(test2)

anno<-c("#DBDB29","#BCBD22","#82853B","#FFC156","#FFB579","#FE870D","#FDA440","#ED8F47","#B15928","#ED4F50","#E4201F","#B294C7","#865FAB","#D4A55B",
        "#B89B74","#99CD91","#69B764","#309343","#A6CEE3","#68A6CD","#569EA4","#2A7FB7","#595959","#BD0A36","#FFB6B0")
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
         gaps_col=c(12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92),
         cellheight=5,
         cellwidth=3,
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.1))),
         annotation_col=df,
         annotation_colors=list(celltype=anno))

# Add clinical information
df<-read.delim('261020_Hannover_sample_table.txt', header=T)
df<-df[,c(1,2,4:11,22,24,27:32)]

combined_list[['disease']]<-df$disease[match(combined_list@meta.data$LuT, df$LuT)]
combined_list[['disease_tissue']]<-df$disease_tissue[match(combined_list@meta.data$LuT, df$LuT)]
combined_list[['classification']]<-df$Classification[match(combined_list@meta.data$LuT, df$LuT)]
combined_list[['other_diagnosis']]<-df$other_diagnosis[match(combined_list@meta.data$LuT, df$LuT)]
combined_list[['site']]<-df$site[match(combined_list@meta.data$LuT, df$LuT)]
combined_list[['immediate']]<-df$Immediate[match(combined_list@meta.data$LuT, df$LuT)]
combined_list[['age']]<-df$Age[match(combined_list@meta.data$LuT, df$LuT)]
combined_list[['sex']]<-df$sex[match(combined_list@meta.data$LuT, df$LuT)]
combined_list[['BMI']]<-df$BMI[match(combined_list@meta.data$LuT, df$LuT)]
combined_list[['blood_group']]<-df$Blood_group[match(combined_list@meta.data$LuT, df$LuT)]
combined_list[['CMV_vaccination']]<-df$IgG_CMV_Empf[match(combined_list@meta.data$LuT, df$LuT)]
combined_list[['FEV1']]<-df$FEV1[match(combined_list@meta.data$LuT, df$LuT)]
combined_list[['FVC']]<-df$FVC1[match(combined_list@meta.data$LuT, df$LuT)]
combined_list[['smoking']]<-df$Nicotine_abuse[match(combined_list@meta.data$LuT, df$LuT)]