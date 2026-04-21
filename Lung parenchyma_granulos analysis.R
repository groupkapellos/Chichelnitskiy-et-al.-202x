# Set directory
setwd('~/hannover/analysis/Theo_analysis')

# All granulocytes
DefaultAssay(object=combined_list)<-"integrated"
Idents(combined_list)<-combined_list$integrated_snn_res.0.4
granulos<-subset(combined_list, idents=c(7,17,16))

Idents(granulos)<-'doublets_individual'
granulos<-subset(granulos, idents='Singlet')

# Find variable genes
granulos<-FindVariableFeatures(object=granulos, selection.method="vst", nfeatures=2000)

# Scale data
granulos<-ScaleData(object=granulos)

# Perform linear reduction    
granulos<-RunPCA(object=granulos, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=granulos, ndims=50)

# Find k-nearest neighbours
DefaultAssay(granulos)<-'integrated'
granulos<-FindNeighbors(object=granulos, reduction='pca', dims=1:20)

# Run UMAP
granulos<-RunUMAP(object=granulos, reduction='pca', dims=1:20)
DimPlot(object=granulos, reduction='umap', group.by='general_final_annotation', label=F, pt.size=0.1)+scale_color_manual(values=c('#FEB70D','#FDA440','#ED8F47','#B15928'))


## Neutrophils
DefaultAssay(object=granulos)<-"integrated"
Idents(granulos)<-granulos$general_final_annotation
neutros<-subset(granulos, idents=unique(granulos$general_final_annotation)[c(2,4)])

# Find variable genes
neutros<-FindVariableFeatures(object=neutros, selection.method="vst", nfeatures=2000)

# Scale data
neutros<-ScaleData(object=neutros)

# Perform linear reduction    
neutros<-RunPCA(object=neutros, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=neutros, ndims=50)

# Find k-nearest neighbours
DefaultAssay(neutros)<-'integrated'
neutros<-FindNeighbors(object=neutros, reduction='pca', dims=1:25)

# Find clusters
neutros<-FindClusters(object=neutros, resolution=0.4, group.singletons=T)

# Run UMAP
neutros<-RunUMAP(object=neutros, reduction='pca', dims=1:25)
DimPlot(object=neutros, reduction='umap', label=F, pt.size=0.1) + scale_color_jco()

# Find cluster markers
Idents(neutros)<-neutros$disease
neutros<-subset(neutros, idents=c("Emphysema", "Fibrosis", "Tu.free"))

DefaultAssay(neutros)<-'RNA'

Idents(neutros)<-neutros$integrated_snn_res.0.4

results.neutros<-FindAllMarkers(object=neutros, min.pct=0.20, log=0.25)
results.neutros<-results.neutros[results.neutros$p_val_adj<=0.05 & results.neutros$avg_log2FC>0,]
results.neutros<-results.neutros %>% arrange(cluster, -avg_log2FC)
results.neutros5<-results.neutros %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)

# Plot DE genes between clusters
test<-neutros@assays$RNA@data  
test<-test[rowSums(as.matrix(test))>0,]

tmp2<-data.frame(rownames(test))
colnames(tmp2)<-'genes'
neutros@meta.data$heatmap<-paste(neutros@meta.data$disease, sep='_', neutros@meta.data$integrated_snn_res.0.4)

for (i in sort(as.character(unique(neutros@meta.data$heatmap)),decreasing=F)) {
  tmp<-test[,colnames(test) %in% rownames(neutros@meta.data[neutros@meta.data$heatmap==i,])]
  tmp<-as.data.frame(tmp)
  tmp<-as.data.frame(rowMeans(tmp))
  colnames(tmp)<-i
  tmp2<-cbind(tmp2,tmp)
}  
tmp2$genes<-NULL

test2<-tmp2[rownames(tmp2) %in% unique(results.neutros5$gene),]
test2<-test2[unique(results.neutros5$gene),]
test2<-test2[,c(11,1,7,
                12,2,8,
                13,3,9,
                14,4,10,
                15,5)]

df<-data.frame(disease=c(rep(c('Tumor-free', 'Emphysema', 'Fibrosis'), 4), 'Tumor-free', 'Emphysema'),
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
         gaps_col=c(3,6,9,12,14),
         annotation_col=df,
         breaks=seq(-2, 2, by=0.05),
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.05))))


## Mast cells
DefaultAssay(object=granulos)<-"integrated"
Idents(granulos)<-granulos$general_final_annotation
mast<-subset(granulos, idents=unique(granulos$general_final_annotation)[1])

# Find variable genes
mast<-FindVariableFeatures(object=mast, selection.method="vst", nfeatures=2000)

# Scale data
mast<-ScaleData(object=mast)

# Perform linear reduction    
mast<-RunPCA(object=mast, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=mast, ndims=50)

# Find k-nearest neighbours
DefaultAssay(mast)<-'integrated'
mast<-FindNeighbors(object=mast, reduction='pca', dims=1:20)

# Find clusters
mast<-FindClusters(object=mast, resolution=0.2, group.singletons=T)

# Run UMAP
mast<-RunUMAP(object=mast, reduction='pca', dims=1:20)
DimPlot(object=mast, reduction='umap', label=F, pt.size=0.1) + scale_color_uchicago()

# Find cluster markers
Idents(mast)<-mast$disease
mast<-subset(mast, idents=c("Emphysema", "Fibrosis", "Tu.free"))

DefaultAssay(mast)<-'RNA'

Idents(mast)<-mast$integrated_snn_res.0.2

results.mast<-FindAllMarkers(object=mast, min.pct=0.20, log=0.25)
results.mast<-results.mast[results.mast$p_val_adj<=0.05 & results.mast$avg_log2FC>0,]
results.mast<-results.mast %>% arrange(cluster, -avg_log2FC)
results.mast5<-results.mast %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)

# Plot DE genes between clusters
test<-mast@assays$RNA@data  

test<-test[rowSums(as.matrix(test))>0,]
tmp2<-data.frame(rownames(test))
colnames(tmp2)<-'genes'
mast@meta.data$heatmap<-paste(mast@meta.data$disease, sep='_', mast@meta.data$integrated_snn_res.0.2)

for (i in sort(as.character(unique(mast@meta.data$heatmap)),decreasing=F)) {
  tmp<-test[,colnames(test) %in% rownames(mast@meta.data[mast@meta.data$heatmap==i,])]
  tmp<-as.data.frame(tmp)
  tmp<-as.data.frame(rowMeans(tmp))
  colnames(tmp)<-i
  tmp2<-cbind(tmp2,tmp)
}  
tmp2$genes<-NULL

test2<-tmp2[rownames(tmp2) %in% unique(results.mast5$gene),]
test2<-test2[unique(results.mast5$gene),]
test2<-test2[,c(7,1,4,
                8,2,5,
                9,3,6)]

df<-data.frame(disease=c(rep(c('Tumor-free', 'Emphysema', 'Fibrosis'), 3)),
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
         gaps_col=c(3,6),
         annotation_col=df,
         breaks=seq(-2, 2, by=0.05),
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.05))))


# AUCell analysis of BALF neutrophils
library(AUCell)
library(readxl)
library(doMC)
library(doRNG)

tmp<-as.matrix(neutros@assays$RNA@data)
tmp2<-rowSums(tmp)[rowSums(tmp)>3]
tmp<-tmp[row.names(tmp) %in% names(tmp2),]
dim(tmp)

cells_rankings<-AUCell_buildRankings(tmp)
cells_rankings

granules<-list(azurophil=c('CD63','CD68','MPO','BPI','DEFA3','DEFA4','ELANE','AZU1','CTSG','PRTN3','NSP4','CTSC','SERPINA1','LYZ'),
               specific=c('ITGAM','ITGB2','CEACAM8','CD177','NOX2','SCAMP1','SERPINA1','LYZ','B2M','MMP1','MMP8','MMP9','MMP13','MMP2','HP','CAMP','LTF','LCN2','PTX3','SLPI','OLFM4'),
               gelatinase=c('ITGAM','ITGB2','CD177','NOX2','MMP25','FPR1','SCAMP1','VAMP2','LYZ','ARG1','MMP2','MMP9','FCN1'),
               secretory=c('MME','ITGAM','ITGB2','FUT4','ALP','FCGR3B','CR1','NOX2','MMP25','SLC11A2','SCAMP1','VAMP2','SERPINA1','HP'))

cells_AUC<-AUCell_calcAUC(granules, cells_rankings, aucMaxRank=ceiling(0.05*nrow(cells_rankings)), normAUC=T)
cells_AUC

res<-as.data.frame(t(getAUC(cells_AUC)))
res$cluster<-neutros$integrated_snn_res.0.4
res<-melt(res, id.vars=c('cluster'))

res %>% group_by(cluster,variable) %>% summarize(Mean=mean(value))->df.tmp
res<-dcast(formula=variable~cluster, data=df.tmp, value.var='Mean')
rownames(res)<-res$variable
res$variable<-NULL

pheatmap(res[,1:5],
         clustering_distance_rows='correlation',
         clustering_distance_cols='correlation',
         scale='row',
         show_rownames=T,
         show_colnames=T,
         cluster_rows=F,
         cluster_cols=F,
         cellheight=30,
         cellwidth=35,
         breaks=seq(0, 1, by=0.05),
         color=colorRampPalette(brewer.pal(n=9, name="Blues"))(length(seq(0, 1, by=0.05))))
