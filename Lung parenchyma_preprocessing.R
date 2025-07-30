# Remove existing variables from R memory
rm(list=ls())

# Load packages
list.of.packages<-c("gdata","data.table","ggplot2","scales","useful","Seurat","readr","dplyr","RColorBrewer","Matrix","ggExtra","backports","jsonlite","ggrepel")
new.packages<-list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)>0) install.packages(new.packages,suppressUpdates=TRUE)

list.of.bioc.packages<-c("genefilter","clusterProfiler","DOSE","org.Hs.eg.db")

new.packages.bioc<-list.of.bioc.packages[!(list.of.bioc.packages %in% installed.packages()[,"Package"])]

source("https://bioconductor.org/biocLite.R")
if(length(new.packages.bioc)>0) biocLite(new.packages.bioc,suppressUpdates=TRUE)
lapply(c(list.of.packages,list.of.bioc.packages,"SingleR"), require, character.only=TRUE)

# Set working directory
dir<-"~/hannover/alignment/2020-05-22/output/results/samples"
setwd(dir)

# Create Seurat objects for each pool   
sample217<-paste('Sample', sep='_', 204:211)
sample222<-paste('Sample', sep='_', 212:215)
sample225<-paste('Sample', sep='_', 216:219)
sample241<-paste('Sample', sep='_', 220:222)
sample247<-paste('Sample', sep='_', 232:233)
sample264<-paste('Sample', sep='_', 234:238)
sample266<-paste('Sample', sep='_', 239:246)
sample290<-paste('Sample', sep='_', 259:266)
sample294<-paste('Sample', sep='_', 267:272)
sample314<-paste('Sample', sep='_', c(275:277, 280:285))     
sample349<-paste('Sample', sep='_', 298:303)
sample383<-paste('Sample', sep='_', 312:317)
sample410<-paste('Sample', sep='_', 318:326)
sample421<-paste('Sample', sep='_', 336:342)
sample437<-paste('Sample', sep='_', 350:354)
sample433<-paste('Sample', sep='_', c(362:365,367,369:370))
sample500<-paste('Sample', sep='_', 377:380)
sample519<-paste('Sample', sep='_', c(396,398,399,401)) 

list_seurat1<-list()
list_seurat1<-foreach(i=c(sample217,sample222,sample225,sample241,sample247,sample264,
                          sample266,sample290,sample294,sample314,sample349,sample383,
                          sample410,sample421,sample433,sample437,sample500,sample519), .packages=c("Matrix","Seurat")) %dopar% {
  mtx<-readMM(paste(dir, i, "umi/matrix.mtx", sep="/"))
  genes<-read.delim(paste(dir, i, "umi/genes.tsv" ,sep="/"), header=F, stringsAsFactors=F)
  bc<-read.delim(paste(dir, i, "umi/barcodes.tsv", sep="/"), header=F, stringsAsFactors=F)
  mtx@Dimnames[[1]]<-as.character(genes$V1)
  mtx@Dimnames[[2]]<-as.character(bc$V1)
  mtx<-mtx[!grepl(c("^MT-RNR1|^MT-RNR2"), mtx@Dimnames[[1]]),]
  seurat<-CreateSeuratObject(counts=mtx, min.cells=5, min.features=50, project=i)
}

names(list_seurat1)<-c(sample217,sample222,sample225,sample241,sample247,sample264,
                       sample266,sample290,sample294,sample314,sample349,sample383,
                       sample410,sample421,sample433,sample437,sample500,sample519)

names(list_seurat1)<-gsub('Sample', 'Pool', names(list_seurat1))

# Merge pools into one Seurat object
keep_sample<-unlist(lapply(list_seurat1, function(x){
  ncol(x)>=200
}))

length(keep_sample[keep_sample==T])

seurat1<-merge(list_seurat1[[1]], y=list_seurat1[2:length(list_seurat1)],
                 add.cell.ids=names(keep_sample))

seurat1$orig.ident<-gsub('Sample', 'Pool', seurat1$orig.ident)
seurat1

# Add alignment data
dir<-"~/hannover/alignment/2020-05-22/output/results/logs/dropseq_tools/"
setwd(dir)

align_stats_list<-list()
for(i in c(sample217,sample222,sample225,sample241,sample247,sample264,
           sample266,sample290,sample294,sample314,sample349,sample383,
           sample410,sample421,sample433,sample437,sample500,sample519)){
  
    rd.stats<-read_delim(paste0(dir, i, '_rna_metrics.txt'), delim='\t', skip=6)
    rd.stats$Cell_ID<-paste0(i, '_',  rd.stats$SAMPLE)
    align_stats_list[[i]]<-rd.stats
}

align_stats_df<-as.data.frame(do.call(bind_rows, align_stats_list))
rownames(align_stats_df)<-align_stats_df$Cell_ID
rownames(align_stats_df)<-gsub('Sample', 'Pool', rownames(align_stats_df))
align_stats_df<-align_stats_df[rownames(align_stats_df) %in% colnames(seurat1@assays$RNA@counts), c(21, 22)]

identical(rownames(align_stats_df), rownames(seurat1@meta.data))
tmp<-align_stats_df
tmp<-merge(seurat1@meta.data, tmp, by=0)
rownames(tmp)<-tmp$Row.names
tmp$Row.names<-NULL
seurat1@meta.data<-tmp

# Set working directory       
dir<-"~/Medimmune_cancer/alignment/2020-08-14_402-403_405-460/output/results/samples"
setwd(dir)

sample426<-paste('sample', sep='_', 454:460)     
sample580<-paste('sample', sep='_', 436:443)

list_seurat2<-list()
list_seurat2<-foreach(i=c(sample426,sample580), .packages=c("Matrix","Seurat")) %dopar% {
                            mtx<-readMM(paste(dir, i, "umi/matrix.mtx", sep="/"))
                            genes<-read.delim(paste(dir, i, "umi/genes.tsv", sep="/"), header=F, stringsAsFactors=F)
                            bc<-read.delim(paste(dir, i, "umi/barcodes.tsv", sep="/"), header=F, stringsAsFactors=F)
                            mtx@Dimnames[[1]]<-as.character(genes$V1)
                            mtx@Dimnames[[2]]<-as.character(bc$V1)
                            mtx<-mx[!grepl(c("^MT-RNR1|^MT-RNR2"), mtx@Dimnames[[1]]),]
                            seurat<-CreateSeuratObject(counts=mtx, min.cells=5, min.features=50, project=i)
}

names(list_seurat2)<-c(sample426, sample580)
names(list_seurat2)<-gsub('sample', 'Pool', names(list_seurat2))

# Merge pools into one Seurat object
keep_sample<-unlist(lapply(list_seurat2, function(x){
  ncol(x)>200
}))

length(keep_sample[keep_sample==T])

seurat2<-merge(list_seurat2[[1]], y=list_seurat2[2:length(list_seurat2)],
                 add.cell.ids=names(keep_sample))

seurat2$orig.ident<-gsub('sample', 'Pool', seurat2$orig.ident)
seurat2

# Add alignment data
dir<-"~/Medimmune_cancer/alignment/2020-08-14_402-403_405-460/output/results/logs/dropseq_tools/"
setwd(dir)

align_stats_list<-list()
for(i in c(sample426, sample580)){
  
  rd.stats<-read_delim(paste0(dir, i, '_rna_metrics.txt'), delim='\t', skip=6)
  rd.stats$Cell_ID<-paste0(i, '_', rd.stats$SAMPLE)
  align_stats_list[[i]]<-rd.stats
}

align_stats_df<-as.data.frame(do.call(bind_rows, align_stats_list))
rownames(align_stats_df)<-align_stats_df$Cell_ID
rownames(align_stats_df)<-gsub('sample', 'Pool', rownames(align_stats_df))
align_stats_df<-align_stats_df[rownames(align_stats_df) %in% colnames(seurat2@assays$RNA@counts), c(21, 22)]

identical(rownames(align_stats_df), rownames(seurat2@meta.data))
tmp<-align_stats_df
tmp<-merge(seurat2@meta.data, tmp, by=0)
rownames(tmp)<-tmp$Row.names
tmp$Row.names<-NULL
seurat2@meta.data<-tmp

# Set working directory
dir<-"~/Medimmune_cancer/alignment/2020-08-16_423-430_477-483/output/results/samples/"
setwd(dir)

sample553<-c(423:426)

list_seurat3<-list()
list_seurat3<-foreach(i=sample553, .packages=c("Matrix","Seurat")) %dopar% {
                            mtx<-readMM(paste(dir, i, "umi/matrix.mtx",sep="/"))
                            genes<-read.delim(paste(dir, i, "umi/genes.tsv", sep="/"), header=F, stringsAsFactors=F)
                            bc<-read.delim(paste(dir, i, "umi/barcodes.tsv", sep="/"), header=F, stringsAsFactors=F)
                            mtx@Dimnames[[1]]<-as.character(genes$V1)
                            mtx@Dimnames[[2]]<-as.character(bc$V1)
                            mtx<-mtx[!grepl(c("^MT-RNR1|^MT-RNR2"), mtx@Dimnames[[1]]),]
                            seurat<-CreateSeuratObject(counts=mtx, min.cells=5, min.features=50, project=paste('Sample', sep='_', i))
                          }

names(list_seurat3)<-paste('Pool', sep='_', c(423:426))

# Merge pools into one Seurat object
keep_sample<-unlist(lapply(list_seurat3, function(x){
  ncol(x)>200
}))

length(keep_sample[keep_sample==T])

seurat3<-merge(list_seurat3[[1]], y=list_seurat3[2:length(list_seurat3)],
                 add.cell.ids=names(keep_sample))

seurat3$orig.ident<-gsub('Sample', 'Pool', seurat3$orig.ident)
seurat3

# Add alignment data
dir<-"~/Medimmune_cancer/alignment/2020-08-16_423-430_477-483/output/results/logs/dropseq_tools/"
setwd(dir)

align_stats_list<-list()
for(i in sample553){
  
  rd.stats<-read_delim(paste0(dir, i, '_rna_metrics.txt'), delim='\t', skip=6)
  rd.stats$Cell_ID<-paste0(i, '_', rd.stats$SAMPLE)
  align_stats_list[[i]]<-rd.stats
}

align_stats_df<-as.data.frame(do.call(bind_rows, align_stats_list))
rownames(align_stats_df)<-align_stats_df$Cell_ID
rownames(align_stats_df)<-paste('Pool', sep='_', rownames(align_stats_df))
align_stats_df<-align_stats_df[rownames(align_stats_df) %in% colnames(seurat3@assays$RNA@counts), c(21, 22)]

identical(rownames(align_stats_df), rownames(seurat3@meta.data))
tmp<-align_stats_df
tmp<-merge(seurat3@meta.data, tmp, by=0)
rownames(tmp)<-tmp$Row.names
tmp$Row.names<-NULL
seurat3@meta.data<-tmp

# Set working directory
dir<-"~/hannover/alignment/2020-10-19/output/results/samples/"
setwd(dir)

sample736<-c(487:490) 

list_seurat4<-list()
list_seurat4<-foreach(i=sample736, .packages=c("Matrix","Seurat")) %dopar% {
  mtx<-readMM(paste(dir, i, "umi/matrix.mtx",sep="/"))
  genes<-read.delim(paste(dir, i, "umi/genes.tsv", sep="/"), header=F, stringsAsFactors=F)
  bc<-read.delim(paste(dir, i, "umi/barcodes.tsv", sep="/"), header=F, stringsAsFactors=F)
  mtx@Dimnames[[1]]<-as.character(genes$V1)
  mtx@Dimnames[[2]]<-as.character(bc$V1)
  mtx<-mtx[!grepl(c("^MT-RNR1|^MT-RNR2"), mtx@Dimnames[[1]]),]
  seurat<-CreateSeuratObject(counts=mtx, min.cells=5, min.features=50, project=paste('Sample', sep='_', i))
}
 
# Merge pools into one Seurat object
keep_sample<-unlist(lapply(list_seurat4, function(x){
  ncol(x)>200
}))

length(keep_sample[keep_sample==T])

seurat4<-merge(list_seurat4[[1]], y=list_seurat4[2:length(list_seurat4)],
               add.cell.ids=names(keep_sample))

seurat4$orig.ident<-gsub('Sample', 'Pool', seurat4$orig.ident)
seurat4

# Add alignment data
dir<-"~/hannover/alignment/2020-10-19/output/results/logs/dropseq_tools/"
setwd(dir)

align_stats_list<-list()
for(i in sample736){
  
  rd.stats<-read_delim(paste0(dir, i, '_rna_metrics.txt'), delim='\t', skip=6)
  rd.stats$Cell_ID<-paste0(i, '_', rd.stats$SAMPLE)
  align_stats_list[[i]]<-rd.stats
}

align_stats_df<-as.data.frame(do.call(bind_rows, align_stats_list))
rownames(align_stats_df)<-align_stats_df$Cell_ID
rownames(align_stats_df)<-paste('Pool', sep='_', rownames(align_stats_df))
align_stats_df<-align_stats_df[rownames(align_stats_df) %in% colnames(seurat4@assays$RNA@counts),  c(21, 22)]

identical(rownames(align_stats_df), rownames(seurat4@meta.data))
tmp<-align_stats_df
tmp<-merge(seurat4@meta.data, tmp, by=0)
rownames(tmp)<-tmp$Row.names
tmp$Row.names<-NULL
seurat4@meta.data<-tmp

# Merge runs
combined<-merge(x=seurat1, y=c(seurat2, seurat3, seurat4))

# Decide on inflection point
df<-data.frame(row.names=rownames(combined@meta.data),
               cell=rownames(combined@meta.data),
               Gene=combined$nFeature_RNA,
               UMI=combined$nCount_RNA)

df<-df[order(df$UMI, decreasing=T),]
df$index<-1:nrow(df)

ggplot(df, aes(Gene, UMI))+
  geom_point(size=0.1, show.legend=F)+
  scale_y_log10()+
  scale_x_log10()+
  geom_abline(slope=1, color='red')+
  geom_vline(aes(xintercept=100), color='orange')+
  geom_vline(aes(xintercept=150), color='brown')+
  ylab('UMI numbers')+
  xlab('Gene numbers')+
  theme_bw()

# Filter on high number of genes
quantile(combined$nFeature_RNA, probs=c(0.25 ,0.5, 0.75, 0.95, 0.99))

combined<-subset(combined, subset=nFeature_RNA>=200)
combined<-subset(combined, subset=nFeature_RNA<=2000)

# Set new directory
dir<-"~/hannover/analysis"
setwd(dir)

# Remove pools with very few cells
sampletable<-read.delim('sampletable.csv', sep=';')

hist(table(combined$orig.ident), 
     main='Distribution of pool size', 
     xlab='Cells (100 cells/bin)', 
     breaks=45,
     col='brown')

numbers<-combined@meta.data %>% group_by(orig.ident) %>% tally()
numbers<-numbers[numbers$n>100,]
numbers$LuT<-sampletable$ID[match(numbers$orig.ident, sampletable$Pool)]

Idents(combined)<-'orig.ident'
combined<-subset(combined, idents=numbers$orig.ident)

# Annotate the samples
combined$LuT<-sampletable$ID[match(combined$orig.ident, sampletable$Pool)]

gc(verbose=T)

# Remove cells due to wrong TSO primer
cell.type.tmp<-row.names(combined@meta.data)
cell.type.tmp<-cell.type.tmp[grep("_ATGGG|_GGG", cell.type.tmp, invert=T)]
cell.type.hamming<-NULL
for (i in 1:length(cell.type.tmp)) {
  tmp.new<-as.matrix(stringDist(c(str_sub(cell.type.tmp[i], -12, -8)), "ATGGG"), method="hamming")[1, 2]
  cell.type.hamming<-c(cell.type.hamming, tmp.new)
}
cell.type.tmp<-cell.type.tmp[cell.type.hamming >= 2]

combined$wrong_barcodes<-'yes'
combined@meta.data[rownames(combined@meta.data) %in% cell.type.tmp,"wrong_barcodes"]<-'no'
combined<-subset(combined, subset=wrong_barcodes=='no')

# Filter on mitochondrial gene content
combined[["percent.mt"]]<-PercentageFeatureSet(object=combined, pattern="^MT-")
combined<-subset(combined, subset=percent.mt<5)

# Filter on intronic read content
summary(combined$PCT_INTRONIC_BASES)
summary(combined$PCT_INTERGENIC_BASES)

combined<-subset(combined, subset=PCT_INTRONIC_BASES<0.3)

# QC analysis
par(mfrow=c(1, 2))
FeatureScatter(object=combined, feature1="nCount_RNA", feature2="percent.mt")
FeatureScatter(object=combined, feature1="nFeature_RNA", feature2="nCount_RNA")

VlnPlot(object=combined, features=c("nFeature_RNA", "nCount_RNA", "percent.mt", "PCT_INTRONIC_BASES"), pt.size=0, log=F, ncol=1)

# Final removal of samples due to low numbers
Idents(combined)<-'LuT'
combined<-subset(combined, idents=c('LuT_247','LuT_264','LuT_266','LuT_290','LuT_294','LuT_314','LuT_349',
                                    'LuT_383','LuT_410','LuT_421','LuT_433','LuT_437','LuT_500','LuT_553',
                                    'LuT_580','LuT_736'))

# Integrate the dataset
data.list<-SplitObject(combined, split.by="LuT")

for (i in 1:length(data.list)) {
  data.list[[i]]<-NormalizeData(data.list[[i]], verbose=FALSE)
  data.list[[i]]<-FindVariableFeatures(data.list[[i]], selection.method="vst", nfeatures=2000, verbose=FALSE)
}

anchors<-FindIntegrationAnchors(object.list=data.list, dims=1:30)
combined_list<-IntegrateData(anchorset=anchors, dims=1:30)

# Scale data
DefaultAssay(object=combined_list)<-"integrated"
combined_list<-ScaleData(object=combined_list)

# Perform linear reduction    
combined_list<-RunPCA(object=combined_list, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=combined_list, ndims=50)

# Find k-nearest neighbours
combined_list<-FindNeighbors(object=combined_list, reduction='pca', dims=1:50)

# Find clusters
DefaultAssay(object=combined_list)<-"integrated"
combined_list<-FindClusters(object=combined_list, resolution=0.6, group.singletons=T)

# Run UMAP
combined_list<-RunUMAP(object=combined_list, reduction='pca', dims=1:50)
DimPlot(object=combined_list, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=combined_list) <- "RNA"
results<-all.markers(object=combined_list, min.pct=0.25, log=0.40)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=avg_logFC)

DotPlot(object=combined_list, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis()

# Remove uninformative cluster (Cluster 0)
Idents(combined_list)<-combined_list$integrated_snn_res.0.6
combined_list<-subset(combined_list, idents=c(1:21))

# Find variable genes
combined_list<-FindVariableFeatures(object=combined_list, selection.method="vst", nfeatures=2000)

# Scale data
combined_list<-ScaleData(object=combined_list)

# Perform linear reduction    
combined_list<-RunPCA(object=combined_list, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=combined_list, ndims=50)

# Find k-nearest neighbours
combined_list<-FindNeighbors(object=combined_list, reduction='pca', dims=1:30)

# Find clusters
combined_list<-FindClusters(object=combined_list, resolution=0.4, group.singletons=T)

# Run UMAP
combined_list<-RunUMAP(object=combined_list, reduction='pca', dims=1:30)
DimPlot(object=combined_list, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=combined_list)<-"RNA"
results<-all.markers(object=combined_list, min.pct=0.20, log=0.40)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=avg_logFC)

DotPlot(object=combined_list, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis()
                        
# Identify doublets(run on every sample separately)
sweep.res.list_combined<-paramSweep_v3(combined, PCs=1:30, sct=FALSE)
sweep.stats_combined<-summarizeSweep(sweep.res.list_combined, GT=FALSE)
bcmvn_combined<-find.pK(sweep.stats_combined)
bcmvn_combined[bcmvn_combined$BCmetric==max(bcmvn_combined$BCmetric),]

homotypic.prop<-modelHomotypic(combined@meta.data$RNA_snn_res.0.4)  
nExp_poi<-round(0.08*length(combined@meta.data$RNA_snn_res.0.4))  
nExp_poi.adj<-round(nExp_poi*(1-homotypic.prop))

combined<-doubletFinder_v3(combined, PCs=1:30, pK=0.07, nExp=nExp_poi, reuse.pANN=FALSE, sct=FALSE)
combined<-doubletFinder_v3(combined, PCs=1:30, pK=0.07, nExp=nExp_poi.adj, reuse.pANN="pANN_0.25_0.07_12492", sct=FALSE)