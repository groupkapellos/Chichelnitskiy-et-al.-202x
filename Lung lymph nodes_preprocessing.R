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
sample329<-paste('Sample', sep='_', c(286:288, 290, 292))
sample358<-paste('Sample', sep='_', 304:309)
sample412<-paste('Sample', sep='_', 327:335)
sample424<-paste('Sample', sep='_', 343:349)
sample436<-paste('Sample', sep='_', 371:376)
sample503<-paste('Sample', sep='_', 381:387)

list_seurat1<-list()
list_seurat1<-foreach(i=c(sample329,sample358,sample412,sample424,sample436,sample503), .packages=c("Matrix","Seurat")) %dopar% {
  mtx<-readMM(paste(dir ,i, "umi/matrix.mtx", sep="/"))
  genes<-read.delim(paste(dir,i,"umi/genes.tsv", sep="/"), header=F, stringsAsFactors=F)
  bc<-read.delim(paste(dir,i,"umi/barcodes.tsv", sep="/"), header=F, stringsAsFactors=F)
  mtx@Dimnames[[1]]<-as.character(genes$V1)
  mtx@Dimnames[[2]]<-as.character(bc$V1)
  mtx<-mtx[!grepl(c("^MT-RNR1|^MT-RNR2"), mtx@Dimnames[[1]]),]
  seurat<-CreateSeuratObject(counts=mtx, min.cells=5, min.features=50, project=i)
}

names(list_seurat1)<-c(sample329,sample358,sample412,sample424,sample436,sample503)

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
for(i in c(sample329,sample358,sample412,sample424,sample436,sample503)){
  
    rd.stats<-read_delim(paste0(dir, i, '_rna_metrics.txt'), delim='\t', skip=6)
    rd.stats$Cell_ID<-paste0(i, '_', rd.stats$SAMPLE)
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
dir<-"~/hannover/alignment/2020-09-29/output/results/samples/"
setwd(dir)

sample636<-paste('sample',sep='_',  c(470, 472)) 

list_seurat2<-list()
list_seurat2<-foreach(i=sample636, .packages=c("Matrix","Seurat")) %dopar% {
  mtx<-readMM(paste(dir, i, "umi/matrix.mtx",sep="/"))
  genes<-read.delim(paste(dir, i, "umi/genes.tsv",sep="/"), header=F, stringsAsFactors=F)
  bc<-read.delim(paste(dir, i, "umi/barcodes.tsv",sep="/"), header=F, stringsAsFactors=F)
  mtx@Dimnames[[1]]<-as.character(genes$V1)
  mtx@Dimnames[[2]]<-as.character(bc$V1)
  mtx<-mtx[!grepl(c("^MT-RNR1|^MT-RNR2"), mtx@Dimnames[[1]]),]
  seurat<-CreateSeuratObject(counts=mtx, min.cells=5, min.features=50, project=i)
}

names(list_seurat2)<-paste('Pool', sep='_', c(470, 472))

# Merge pools into one Seurat object
keep_sample<-unlist(lapply(list_seurat2, function(x){
  ncol(x)>200
}))

length(keep_sample[keep_sample==T])

seurat2<-merge(list_seurat2[[1]], y=list_seurat2[[2]],
               add.cell.ids=names(keep_sample))

seurat2$orig.ident<-gsub('sample', 'Pool', seurat2$orig.ident)
seurat2

# Add alignment data
dir<-"~/hannover/alignment/2020-09-29/output/results/logs/dropseq_tools/"
setwd(dir)

align_stats_list<-list()
for(i in sample636){
  
  rd.stats<-read_delim(paste0(dir, i, '_rna_metrics.txt'), delim='\t', skip=6)
  rd.stats$Cell_ID<-paste0(i, '_',rd.stats$SAMPLE)
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
dir<-"~/hannover/alignment/2020-10-19/output/results/samples/"
setwd(dir)

sample737<-c(491:494, 496) 

list_seurat3<-list()
list_seurat3<-foreach(i=sample737, .packages=c("Matrix","Seurat")) %dopar% {
  mtx<-readMM(paste(dir ,i, "umi/matrix.mtx",sep="/"))
  genes<-read.delim(paste(dir, i, "umi/genes.tsv",sep="/"), header=F, stringsAsFactors=F)
  bc<-read.delim(paste(dir, i, "umi/barcodes.tsv",sep="/"), header=F, stringsAsFactors=F)
  mtx@Dimnames[[1]]<-as.character(genes$V1)
  mtx@Dimnames[[2]]<-as.character(bc$V1)
  mtx<-mtx[!grepl(c("^MT-RNR1|^MT-RNR2"), mtx@Dimnames[[1]]),]
  seurat<-CreateSeuratObject(counts=mtx, min.cells=5, min.features=50, project=paste('Sample', sep='_', i))
}

names(list_seurat3)<-paste('Pool', sep='_', c(491:494, 496))

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
dir<-"~/hannover/alignment/2020-10-19/output/results/logs/dropseq_tools/"
setwd(dir)

align_stats_list<-list()
for(i in sample737){
  
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

# Merge runs
combined<-merge(x=seurat1, y=c(seurat2, seurat3))

# Decide on inflection point
df<-data.frame(row.names=rownames(combined@meta.data),
               cell=rownames(combined@meta.data),
               Gene=combined$nFeature_RNA,
               UMI=combined$nCount_RNA)

df<-df[order(df$UMI, decreasing=T),]
df$index<-1:nrow(df)

ggplot(df, aes(index, UMI))+
  geom_point(size=0.1, show.legend=F)+
  scale_y_log10()+
  scale_x_log10()+
  geom_hline(aes(yintercept=150))+
  geom_hline(aes(yintercept=100), color='red')+
  geom_hline(aes(yintercept=50), color='blue')+
  ylab('UMI numbers')+
  xlab('cells sorted by transcript count')+
  theme_bw()

# Filter on high number of genes
quantile(combined$nFeature_RNA, probs=c(0.25, 0.5, 0.75, 0.95, 0.99))

combined<-subset(combined, subset=nFeature_RNA>=200)
#combined<-subset(combined, subset=nFeature_RNA<=2000)

# Set new directory
dir<-"~/hannover/analysis"
setwd(dir)

# Remove pools with very few cells
sampletable<-read.delim('sampletable2.csv', sep=';')

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
  tmp.new<-as.matrix(stringdist(str_sub(cell.type.tmp[i], -12, -8), "ATGGG", method="hamming"))[1]
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

# Integrate the dataset
Idents(combined)<-'LuT'
data.list<-SplitObject(combined, split.by="LuT")

for (i in 1:length(data.list)) {
  data.list[[i]]<-NormalizeData(data.list[[i]], verbose=FALSE)
  data.list[[i]]<-FindVariableFeatures(data.list[[i]], selection.method="vst", nfeatures=2000, verbose=FALSE)
}

anchors<-FindIntegrationAnchors(object.list=data.list, dims=1:30)
integrated.lymph<-IntegrateData(anchorset=anchors, dims=1:30)

# Scale data
DefaultAssay(object=integrated.lymph)<-"integrated"
integrated.lymph<-ScaleData(object=integrated.lymph)

# Perform linear reduction    
integrated.lymph<-RunPCA(object=integrated.lymph, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=integrated.lymph, ndims=50)

# Find k-nearest neighbours
integrated.lymph<-FindNeighbors(object=integrated.lymph, reduction='pca', dims=1:50)

# Find clusters
integrated.lymph<-FindClusters(object=integrated.lymph, resolution=0.6)

# Run UMAP
integrated.lymph<-RunUMAP(object=integrated.lymph, reduction='pca', dims=1:50)
DimPlot(object=integrated.lymph, reduction='umap', label=T, pt.size=0.1)

# Remove clusters with low numbers
Idents(integrated.lymph)<-integrated.lymph$integrated_snn_res.0.6
integrated.lymph<-subset(integrated.lymph, idents=c(0:13))

# Find variable genes
integrated.lymph<-FindVariableFeatures(object=integrated.lymph, selection.method="vst", nfeatures=2000)

# Scale data
integrated.lymph<-ScaleData(object=integrated.lymph)

# Perform linear reduction    
integrated.lymph<-RunPCA(object=integrated.lymph, npcs=50, verbose=F)

# Decide on which PCs to keep
ElbowPlot(object=integrated.lymph, ndims=50)

# Find k-nearest neighbours
integrated.lymph<-FindNeighbors(object=integrated.lymph, reduction='pca', dims=1:15)

# Find clusters
integrated.lymph<-FindClusters(object=integrated.lymph, resolution=0.3, group.singletons=T)

# Run UMAP
integrated.lymph<-RunUMAP(object=integrated.lymph, reduction='pca', dims=1:15)
DimPlot(object=integrated.lymph, reduction='umap', label=T, pt.size=0.1)

# Find cluster markers
DefaultAssay(object=integrated.lymph)<-"RNA"
results<-all.markers(object=integrated.lymph, min.pct=0.20, log=0.40)
results<-results[results$p_val_adj<=0.05 & results$avg_logFC>0,]
results5<-results %>% group_by(cluster) %>% top_n(n=5, wt=-avg_logFC)

DotPlot(object=integrated.lymph, features=rev(unique(results5$gene)), cols='RdBu', dot.scale=8, dot.min=0.1) + RotatedAxis()