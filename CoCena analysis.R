# Load packages
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(grid)

# Load data
modules<-read.delim('module_genes.txt')
modules<-modules[modules$color %in% c('turquoise','darkorange','lightblue','steelblue','orchid','lightgreen'),]

list1<-list()
for(i in unique(modules$color)){
  tmp<-modules[modules$color==i,]$gene
  list1[[i]]<-tmp
}

## COPD
# Load the datasets
copd<-read.delim('GSE37768.txt')
rownames(copd)<-copd$Symbol
copd<-copd[list1$turquoise,]
copd$Symbol<-NULL
copd<-as.matrix(copd)

# Fibrosis
# Load the datasets
fibrosis<-read.delim('GSE48149.txt')
rownames(fibrosis)<-fibrosis$Symbol
fibrosis<-fibrosis[list1$steelblue,]
fibrosis$Symbol<-NULL
fibrosis<-as.matrix(fibrosis)

# Plot results
sample_info<-read.csv('GSE37768_metadata.csv')
df<-data.frame(disease=sample_info$Stage)
rownames(df)<-colnames(copd)

selected_rows<-c('NDUFB1', 'RAB5A', 'SCARB2', 'TSN', 'MBNL1', 'PINK1', 'PSMC5', 'ALKBH5', 'AKT1', 'NUP98', 'AAAS', 'CDK9', 'SAFB', 'LIMD1', 'LAPTM4B', 'RPS7', 'PSMD10', 'RIPK1', 'NUDT21', 'YWHAZ', 'FXR1', 'EIF4B',
                 'POLR2J', 'YBX1', 'PSME3', 'ADPRH')

rownames(copd)[!rownames(copd) %in% selected_rows]<-""

pheatmap(copd,
         clustering_distance_rows='correlation',
         clustering_distance_cols='correlation',
         scale='row',
         show_colnames=F,
         show_rownames=T,
         cluster_rows=T,
         cluster_cols=F,
         cellheight=2,
         cellwidth=6,
         gaps_col=c(9,19),
         annotation_col=df,
         breaks=seq(-2, 2, by=0.1),
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.1))))

sample_info<-read.csv('GSE48149_metadata.csv')
df<-data.frame(disease=sample_info$X.Sample_source_name_ch1)
rownames(df)<-colnames(fibrosis)
df$disease<-gsub('SSc-PF', 'SSc PF', df$disease)

selected_rows<-c('MRC1', 'C15orf48', 'PLA2G7', 'C1QA', 'C1QC', 'C1QB', 'SLC31A2', 'MARCO', 'H2AFY', 'HEXB', 'SLC16A10', 'LIPA', 'TREM2', 'GM2A', 'ACP5', 
                 'APOC1', 'CYBB', 'EMILIN2', 'LGMN', 'CTSZ', 'MCEMP1', 'OLR1')

selected_rows<-c('KRT13', 'KRT6A', 'KRT6B', 'DSG3', 'SPRR1B', 'SERPINB13', 'KRT15', 'COL7A1', 'KRT5', 'GRHL1')

selected_rows<-c('ZFP36', 'TRIB1', 'JUND', 'NLRP3', 'CCL2', 'IL6', 'CX3CL1', 'RIPK2', 'TICAM1', 'IL1B', 'CSF3', 'TNFAIP3', 'CXCL2', 'PELI1', 'FOS', 'VCAM1',
                 'IRAK2', 'PDE4B', 'NFKB2', 'NFKB1')

selected_rows<-c('SFRP2', 'FKBP10', 'LOXL3', 'MMP11', 'AEBP1', 'COMP', 'LOX', 'FMOD', 'COL5A1', 'COL1A1', 'COL1A2', 'COL3A1', 'COL5A2', 'LOXL1', 'SERPINH1',
                 'EMILIN1')

rownames(fibrosis)[!rownames(fibrosis) %in% selected_rows]<-""

pheatmap(fibrosis[,c(14:22,1:13,31:43,23:30,44:53)],
         clustering_distance_rows = 'correlation',
         clustering_distance_cols = 'correlation',
         scale='row',
         show_colnames=F,
         show_rownames=T,
         cluster_rows=T,
         cluster_cols=F,
         cellheight=1,
         cellwidth=6,
         gaps_col=c(9,22,35,43),
         annotation_col=df,
         breaks=seq(-2, 2, by=0.1),
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.1))))

# Test individual genes
copd<-read.delim('GSE37768.txt')
rownames(copd)<-copd$Symbol
genes<-c('RPL3', 'RPL10', 'RPL11', 'RPL13A', 'RPS3', 'RPS9', 'RPS14', 'RPS19', 'RPS23', 'RPS27', 'CNOT6L', 'MATR3', 'NCL', 'PABPC1', 'HMGB2', 'PTPRC', 'BTG1', 'FYN', 'TSPYL2', 'PRRC2C', 'DDX3X', 'GOLGA4', 'S100A8', 'S100A9', 'S100A12')
copd<-copd[genes,]
copd<-copd[!grepl('NA', rownames(copd)),]
copd$Symbol<-NULL
sample_info<-read.csv('GSE37768_metadata.csv')
df<-data.frame(disease=sample_info$Stage)
rownames(df)<-colnames(copd)

pheatmap(as.matrix(copd),
         clustering_distance_rows = 'correlation',
         clustering_distance_cols = 'correlation',
         scale='row',
         show_colnames=F,
         show_rownames=T,
         cluster_rows=T,
         cluster_cols=F,
         cellheight=10,
         cellwidth=4,
         gaps_col=c(9,19),
         annotation_col=df,
         breaks=seq(-2, 2, by=0.1),
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.1))))
##
fibrosis<-read.delim('GSE48149.txt')
rownames(fibrosis)<-fibrosis$Symbol
genes<-c('SPP1', 'FN1', 'GPNMB', 'CAPG', 'FABP5', 'TREM2', 'CD63', 'ALCAM', 'FABP3', 'GM2A', 'ACSL1', 'ACSL3', 'ACSL4', 'CD36', 'NPC2', 'SOAT1', 'HEXB', 'PSAP', 'DBI', 'LAPTM5', 'LAMP1', 'SCARB2', 'AP2S1', 'AP2M1,
         APOE', 'ANXA2', 'ANXA5', 'TIMP1', 'PFN1', 'ARPC5', 'ACTB', 'ACTG1', 'GAPDH', 'CFL1', 'FMN1', 'SWAP70', 'CYFIP1', 'CAPZB', 'CAPZA1', 'TFRC', 'FTL', 'CTSB', 'CTSD', 'CTSL', 'HSP90B1', 'HSPA5', 'FKBP1A')
fibrosis<-fibrosis[genes,]
fibrosis<-fibrosis[!grepl('NA', rownames(fibrosis)),]
fibrosis$Symbol<-NULL

sample_info<-read.csv('GSE48149_metadata.csv')
df<-data.frame(disease=sample_info$X.Sample_source_name_ch1)
rownames(df)<-colnames(fibrosis)
df$disease<-gsub('SSc-PF', 'SSc PF', df$disease)

pheatmap(as.matrix(fibrosis)[,c(14:22,31:53,23:30,1:13)],
         clustering_distance_rows = 'correlation',
         clustering_distance_cols = 'correlation',
         scale='row',
         show_colnames=F,
         show_rownames=T,
         cluster_rows=T,
         cluster_cols=F,
         cellheight=10,
         cellwidth=4,
         gaps_col=c(9,22,32,40),
         annotation_col=df,
         breaks=seq(-2, 2, by=0.1),
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.1))))