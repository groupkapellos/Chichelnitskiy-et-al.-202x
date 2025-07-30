# Load packages
library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(reshape2)
library(ggplot2)
library(car)

#present in emphysema & fibrosis (>3 cells)
#present in parenchyma & ln (>3 cells)

lung macs (0-5) 0-3 # 5 only in fibrosis, 4 no cycling macs in lns
lung monos (0-5) 0-5
lung cdcs (0-4) 0,1,3 # 4 no cycling dcs in lns
lung pdcs (0-4) 2
lung T (0-10) 0-3, 5-6 # 4 only in emphysema, 8-10 only in emphysema/tumor-free
lung cycling T (0-10) 7
lung B (0-4) 1 
lung plasma (0-4) 2-3 # 0 contamination, 4 only in fibrosis
lung NK (0-3) 0-3

ln macs (0-4) 0-1 # 2-4 only in emphysema
ln monos (0-2) 0 # 1-2 only in emphysema
ln cdcs (0-2) 0 # 2 only in fibrosis/contamination
ln pdcs (0-2) 1
ln T (0-6) 0-5
ln cycling T (0-6) 6
ln B (0-8) 1
ln plasma (0-8) 0, 2 # 5-8 only in fibrosis, 3, 4 only in LuT#503
ln NK (0-2) 1 # 0 contamination, 2 only in fibrosis

# Macrophages
Idents(macs)<-'integrated_snn_res.0.3'
macs.pare<-subset(macs, idents=0:3)

Idents(macs.pare)<-'disease'
macs.pare<-subset(macs.pare, idents=c('Emphysema', 'Fibrosis'))
macs.pare$tissue<-'Parenchyma'

Idents(macs)<-'integrated_snn_res.0.3'
macs.ln<-subset(macs, idents=0:1)
macs.ln$tissue<-'Lymphnode'

merge.macs<-merge(macs.pare, macs.ln, add.cell.ids=c('pare', 'ln'))
merge.macs$annotation<-'macs'

Idents(merge.macs)<-'tissue'
genes<-FindMarkers(merge.macs, ident.1='Parenchyma', ident.2='Lymphnode', min.pct=0.20, logfc.threshold=0.25, test.use='MAST', only.pos=F)
genes<-genes[genes$p_val_adj<=0.05,]
genes<-genes %>% arrange(-avg_log2FC)
genes$symbol<-rownames(genes)

GO<-enrichGO(gene=genes[genes$avg_log2FC>0,]$symbol,
             OrgDb=org.Hs.eg.db,
             ont='BP',
             keyType='SYMBOL',
             minGSSize=10,
             maxGSSize=500,
             pvalueCutoff=0.05,
             readable=FALSE)

GO<-simplify(GO, cutoff=0.7, by="p.adjust", select_fun=min)

df1<-as.data.frame(GO)
df1<-df1[df1$Count>=3,]
df1$tissue<-'Parenchyma'

write.csv(df1, 'macs.csv')

# Monocytes
Idents(monos)<-'integrated_snn_res.0.4'
monos.pare<-subset(monos, idents=0:5)

Idents(monos.pare)<-'disease'
monos.pare<-subset(monos.pare, idents=c('Emphysema', 'Fibrosis'))
monos.pare$tissue<-'Parenchyma'

Idents(monos)<-'integrated_snn_res.0.3'
monos.ln<-subset(monos, idents=0)
monos.ln$tissue<-'Lymphnode'

merge.monos<-merge(monos.pare, monos.ln, add.cell.ids=c('pare', 'ln'))
merge.monos$annotation<-'monos'

Idents(merge.monos)<-'tissue'
genes<-FindMarkers(merge.monos, ident.1='Parenchyma', ident.2='Lymphnode', min.pct=0.20, logfc.threshold=0.25, test.use='MAST', only.pos=F)
genes<-genes[genes$p_val_adj<=0.05,]
genes<-genes %>% arrange(-avg_log2FC)
genes$symbol<-rownames(genes)

GO<-enrichGO(gene=genes[genes$avg_log2FC>0,]$symbol,
             OrgDb=org.Hs.eg.db,
             ont='BP',
             keyType='SYMBOL',
             minGSSize=10,
             maxGSSize=500,
             pvalueCutoff=0.05,
             readable=FALSE)

GO<-simplify(GO, cutoff=0.7, by="p.adjust", select_fun=min)

df1<-as.data.frame(GO)
df1<-df1[df1$Count>=3,]
df1$tissue<-'Parenchyma'

GO<-enrichGO(gene=genes[genes$avg_log2FC<0,]$symbol,
             OrgDb=org.Hs.eg.db,
             ont='BP',
             keyType='SYMBOL',
             minGSSize=10,
             maxGSSize=500,
             pvalueCutoff=0.05,
             readable=FALSE)

GO<-simplify(GO, cutoff=0.7, by="p.adjust", select_fun=min)

df2<-as.data.frame(GO)
df2<-df2[df2$Count>=3,]
df2$tissue<-'Lymphnode'

monos<-rbind(df1, df2)
write.csv(monos, 'monos.csv')

# cDCs
Idents(dcs)<-'integrated_snn_res.0.2'
cdcs.pare<-subset(dcs, idents=c(0,1,3))

Idents(cdcs.pare)<-'disease'
cdcs.pare<-subset(cdcs.pare, idents=c('Emphysema', 'Fibrosis'))
cdcs.pare$tissue<-'Parenchyma'

Idents(dcs)<-'integrated_snn_res.0.5'
cdcs.ln<-subset(dcs, idents=0)
cdcs.ln$tissue<-'Lymphnode'

merge.cdcs<-merge(cdcs.pare, cdcs.ln, add.cell.ids=c('pare', 'ln'))
merge.cdcs$annotation<-'cdcs'

Idents(merge.cdcs)<-'tissue'
genes<-FindMarkers(merge.cdcs, ident.1='Parenchyma', ident.2='Lymphnode', min.pct=0.20, logfc.threshold=0.25, test.use='MAST', only.pos=F)
genes<-genes[genes$p_val_adj<=0.05,]
genes<-genes %>% arrange(-avg_log2FC)
genes$symbol<-rownames(genes)

GO<-enrichGO(gene=genes[genes$avg_log2FC>0,]$symbol,
             OrgDb=org.Hs.eg.db,
             ont='BP',
             keyType='SYMBOL',
             minGSSize=10,
             maxGSSize=500,
             pvalueCutoff=0.05,
             readable=FALSE)

GO<-simplify(GO, cutoff=0.7, by="p.adjust", select_fun=min)

df1<-as.data.frame(GO)
df1<-df1[df1$Count>=3,]
df1$tissue<-'Parenchyma'

GO<-enrichGO(gene=genes[genes$avg_log2FC<0,]$symbol,
             OrgDb=org.Hs.eg.db,
             ont='BP',
             keyType='SYMBOL',
             minGSSize=10,
             maxGSSize=500,
             pvalueCutoff=0.05,
             readable=FALSE)

GO<-simplify(GO, cutoff=0.7, by="p.adjust", select_fun=min)

df2<-as.data.frame(GO)
df2<-df2[df2$Count>=3,]
df2$tissue<-'Lymphnode'

cdcs<-rbind(df1, df2)
write.csv(cdcs, 'cdcs.csv')

# pDCs
Idents(dcs)<-'integrated_snn_res.0.2'
pdcs.pare<-subset(dcs, idents=2)

Idents(pdcs.pare)<-'disease'
pdcs.pare<-subset(pdcs.pare, idents=c('Emphysema', 'Fibrosis'))
pdcs.pare$tissue<-'Parenchyma'

Idents(dcs)<-'integrated_snn_res.0.5'
pdcs.ln<-subset(dcs, idents=1)
pdcs.ln$tissue<-'Lymphnode'

merge.pdcs<-merge(pdcs.pare, pdcs.ln, add.cell.ids=c('pare', 'ln'))
merge.pdcs$annotation<-'pdcs'

# T cells
Idents(tcells)<-'integrated_snn_res.0.7'
tcells.pare<-subset(tcells, idents=c(0:3, 5:6))

Idents(tcells.pare)<-'disease'
tcells.pare<-subset(tcells.pare, idents=c('Emphysema', 'Fibrosis'))
tcells.pare$tissue<-'Parenchyma'

Idents(tcells)<-'integrated_snn_res.0.4'
tcells.ln<-subset(tcells, idents=0:5)
tcells.ln$tissue<-'Lymphnode'

merge.tcells<-merge(tcells.pare, tcells.ln, add.cell.ids=c('pare', 'ln'))
merge.tcells$annotation<-'T cells'

Idents(merge.tcells)<-'tissue'
genes<-FindMarkers(merge.tcells, ident.1='Parenchyma', ident.2='Lymphnode', min.pct=0.20, logfc.threshold=0.25, test.use='MAST', only.pos=F)
genes<-genes[genes$p_val_adj<=0.05,]
genes<-genes %>% arrange(-avg_log2FC)
genes$symbol<-rownames(genes)

GO<-enrichGO(gene=genes[genes$avg_log2FC>0,]$symbol,
             OrgDb=org.Hs.eg.db,
             ont='BP',
             keyType='SYMBOL',
             minGSSize=10,
             maxGSSize=500,
             pvalueCutoff=0.05,
             readable=FALSE)

GO<-simplify(GO, cutoff=0.7, by="p.adjust", select_fun=min)

df1<-as.data.frame(GO)
df1<-df1[df1$Count>=3,]
df1$tissue<-'Parenchyma'

GO<-enrichGO(gene=genes[genes$avg_log2FC<0,]$symbol,
             OrgDb=org.Hs.eg.db,
             ont='BP',
             keyType='SYMBOL',
             minGSSize=10,
             maxGSSize=500,
             pvalueCutoff=0.05,
             readable=FALSE)

GO<-simplify(GO, cutoff=0.7, by="p.adjust", select_fun=min)

df2<-as.data.frame(GO)
df2<-df2[df2$Count>=3,]
df2$tissue<-'Lymphnode'

tcells<-rbind(df1, df2)
write.csv(tcells, 'tcells.csv')

# Cycling T cells
Idents(tcells)<-'integrated_snn_res.0.7'
cycltcells.pare<-subset(tcells, idents=7)

Idents(cycltcells.pare)<-'disease'
cycltcells.pare<-subset(cycltcells.pare, idents=c('Emphysema', 'Fibrosis'))
cycltcells.pare$tissue<-'Parenchyma'

Idents(tcells)<-'integrated_snn_res.0.4'
cycltcells.ln<-subset(tcells, idents=6)
cycltcells.ln$tissue<-'Lymphnode'

merge.cycltcells<-merge(cycltcells.pare, cycltcells.ln, add.cell.ids=c('pare', 'ln'))
merge.cycltcells$annotation<-'Cycling T cells'

Idents(merge.cycltcells)<-'tissue'
genes<-FindMarkers(merge.cycltcells, ident.1='Parenchyma', ident.2='Lymphnode', min.pct=0.20, logfc.threshold=0.25, test.use='MAST', only.pos=F)
genes<-genes[genes$p_val_adj<=0.05,]
genes<-genes %>% arrange(-avg_log2FC)
genes$symbol<-rownames(genes)

GO<-enrichGO(gene=genes[genes$avg_log2FC>0,]$symbol,
             OrgDb=org.Hs.eg.db,
             ont='BP',
             keyType='SYMBOL',
             minGSSize=10,
             maxGSSize=500,
             pvalueCutoff=0.05,
             readable=FALSE)

GO<-simplify(GO, cutoff=0.7, by="p.adjust", select_fun=min)

df1<-as.data.frame(GO)
df1<-df1[df1$Count>=3,]
df1$tissue<-'Parenchyma'

GO<-enrichGO(gene=genes[genes$avg_log2FC<0,]$symbol,
             OrgDb=org.Hs.eg.db,
             ont='BP',
             keyType='SYMBOL',
             minGSSize=10,
             maxGSSize=500,
             pvalueCutoff=0.05,
             readable=FALSE)

GO<-simplify(GO, cutoff=0.7, by="p.adjust", select_fun=min)

df2<-as.data.frame(GO)
df2<-df2[df2$Count>=3,]
df2$tissue<-'Lymphnode'

cycltcells<-df1
write.csv(cycltcells, 'cyclingtcells.csv')

# B cells
Idents(bcells)<-'integrated_snn_res.0.2'
bcells.pare<-subset(bcells, idents=1)
bcells.pare$tissue<-'Parenchyma'

Idents(bcells.pare)<-'disease'
bcells.pare<-subset(bcells.pare, idents=c('Emphysema', 'Fibrosis'))

Idents(bcells)<-'integrated_snn_res.0.2'
bcells.ln<-subset(bcells, idents=1)
bcells.ln$tissue<-'Lymphnode'

merge.bcells<-merge(bcells.pare, bcells.ln, add.cell.ids=c('pare', 'ln'))
merge.bcells$annotation<-'bcells'

Idents(merge.bcells)<-'tissue'
genes<-FindMarkers(merge.bcells, ident.1='Parenchyma', ident.2='Lymphnode', min.pct=0.20, logfc.threshold=0.25, test.use='MAST', only.pos=F)
genes<-genes[genes$p_val_adj<=0.05,]
genes<-genes %>% arrange(-avg_log2FC)
genes$symbol<-rownames(genes)

GO<-enrichGO(gene=genes[genes$avg_log2FC>0,]$symbol,
             OrgDb=org.Hs.eg.db,
             ont='BP',
             keyType='SYMBOL',
             minGSSize=10,
             maxGSSize=500,
             pvalueCutoff=0.05,
             readable=FALSE)

GO<-simplify(GO, cutoff=0.7, by="p.adjust", select_fun=min)

df1<-as.data.frame(GO)
df1<-df1[df1$Count>=3,]
df1$tissue<-'Parenchyma'

GO<-enrichGO(gene=genes[genes$avg_log2FC<0,]$symbol,
             OrgDb=org.Hs.eg.db,
             ont='BP',
             keyType='SYMBOL',
             minGSSize=10,
             maxGSSize=500,
             pvalueCutoff=0.05,
             readable=FALSE)

GO<-simplify(GO, cutoff=0.7, by="p.adjust", select_fun=min)

df2<-as.data.frame(GO)
df2<-df2[df2$Count>=3,]
df2$tissue<-'Lymphnode'

bcells<-rbind(df1, df2)
write.csv(bcells, 'bcells.csv')

# Plasma cells
Idents(bcells)<-'integrated_snn_res.0.2'
plasma.pare<-subset(bcells, idents=2:3)

Idents(plasma.pare)<-'disease'
plasma.pare<-subset(plasma.pare, idents=c('Emphysema', 'Fibrosis'))
plasma.pare$tissue<-'Parenchyma'

Idents(bcells)<-'integrated_snn_res.0.2'
plasma.ln<-subset(bcells, idents=c(0,2))
plasma.ln$tissue<-'Lymphnode'

merge.plasma<-merge(plasma.pare, plasma.ln, add.cell.ids=c('pare', 'ln'))
merge.plasma$annotation<-'plasma'

Idents(merge.plasma)<-'tissue'
genes<-FindMarkers(merge.plasma, ident.1='Parenchyma', ident.2='Lymphnode', min.pct=0.20, logfc.threshold=0.25, test.use='MAST', only.pos=F)
genes<-genes[genes$p_val_adj<=0.05,]
genes<-genes %>% arrange(-avg_log2FC)
genes$symbol<-rownames(genes)

GO<-enrichGO(gene=genes[genes$avg_log2FC>0,]$symbol,
             OrgDb=org.Hs.eg.db,
             ont='BP',
             keyType='SYMBOL',
             minGSSize=10,
             maxGSSize=500,
             pvalueCutoff=0.05,
             readable=FALSE)

GO<-simplify(GO, cutoff=0.7, by="p.adjust", select_fun=min)

df1<-as.data.frame(GO)
df1<-df1[df1$Count>=3,]
df1$tissue<-'Parenchyma'

GO<-enrichGO(gene=genes[genes$avg_log2FC<0,]$symbol,
             OrgDb=org.Hs.eg.db,
             ont='BP',
             keyType='SYMBOL',
             minGSSize=10,
             maxGSSize=500,
             pvalueCutoff=0.05,
             readable=FALSE)

GO<-simplify(GO, cutoff=0.7, by="p.adjust", select_fun=min)

df2<-as.data.frame(GO)
df2<-df2[df2$Count>=3,]
df2$tissue<-'Lymphnode'

plasma<-rbind(df1, df2)
write.csv(plasma, 'plasma.csv')

# NK cells
Idents(nks)<-'integrated_snn_res.0.3'
nks.pare<-subset(nks, idents=0:3)

Idents(nks.pare)<-'disease'
nks.pare<-subset(nks.pare, idents=c('Emphysema', 'Fibrosis'))
nks.pare$tissue<-'Parenchyma'

Idents(nks)<-'integrated_snn_res.0.3'
nks.ln<-subset(nks, idents=1)
nks.ln$tissue<-'Lymphnode'

merge.nks<-merge(nks.pare, nks.ln, add.cell.ids=c('pare', 'ln'))
merge.nks$annotation<-'nks'

Idents(merge.nks)<-'tissue'
genes<-FindMarkers(merge.nks, ident.1='Parenchyma', ident.2='Lymphnode', min.pct=0.20, logfc.threshold=0.25, test.use='MAST', only.pos=F)
genes<-genes[genes$p_val_adj<=0.05,]
genes<-genes %>% arrange(-avg_log2FC)
genes$symbol<-rownames(genes)

GO<-enrichGO(gene=genes[genes$avg_log2FC>0,]$symbol,
             OrgDb=org.Hs.eg.db,
             ont='BP',
             keyType='SYMBOL',
             minGSSize=10,
             maxGSSize=500,
             pvalueCutoff=0.05,
             readable=FALSE)

GO<-simplify(GO, cutoff=0.7, by="p.adjust", select_fun=min)

df1<-as.data.frame(GO)
df1<-df1[df1$Count>=3,]
df1$tissue<-'Parenchyma'

GO<-enrichGO(gene=genes[genes$avg_log2FC<0,]$symbol,
             OrgDb=org.Hs.eg.db,
             ont='BP',
             keyType='SYMBOL',
             minGSSize=10,
             maxGSSize=500,
             pvalueCutoff=0.05,
             readable=FALSE)

GO<-simplify(GO, cutoff=0.7, by="p.adjust", select_fun=min)

df2<-as.data.frame(GO)
df2<-df2[df2$Count>=3,]
df2$tissue<-'Lymphnode'

nks<-rbind(df1, df2)
write.csv(nks, 'nks.csv')

# Calculate disease proportions
merge.all<-merge(merge(merge(merge(merge(merge(merge(merge(merge.macs, merge.monos, add.cell.ids=c('macs', 'monos')), merge.cdcs, add.cell.ids=c('macs/monos', 'cdcs')), merge.pdcs, add.cell.ids=c('macs/monos/cdcs', 'pdcs')), merge.tcells, add.cell.ids=c('macs/monos/cdcs/pdcs', 'tcells')), merge.cycltcells, add.cell.ids=c('macs/monos/cdcs/pdcs/tcells', 'cycltcells')), merge.bcells, add.cell.ids=c('macs/monos/cdcs/pdcs/tcells/cycltcells', 'bcells')), merge.plasma, add.cell.ids=c('macs/monos/cdcs/pdcs/tcells/cycltcells/bcells', 'plasma')), merge.nks, add.cell.ids=c('macs/monos/cdcs/pdcs/tcells/cycltcells/bcells/plasma', 'nks'))
table(merge.all$tissue)
table(merge.all$annotation, merge.all$tissue)

Idents(merge.all)<-'tissue'
tmp<-subset(merge.all, idents='Parenchyma')

set.seed(12345)
data<-round(prop.table(x=table(tmp$LuT, tmp$annotation), margin=1), digits=3)

data<-melt(data)
data$disease<-tmp@meta.data$disease[match(data$Var1, tmp@meta.data$LuT)]
data$disease<-factor(data$disease, levels=c('Emphysema','Fibrosis'))
data$Var2<-factor(data$Var2, levels=c('macs', 'monos', 'cdcs', 'pdcs', 'T cells', 'Cycling T cells', 'bcells', 'plasma', 'nks'))

ggplot(data, aes(x=Var2 ,y=value, fill=disease)) +
  stat_boxplot(geom='errorbar', linetype=1) +
  geom_boxplot(outlier.shape=1) +
  ylab('Percentage') +
  theme(text=element_text(size=14),
        axis.title.y=element_text(size=14,face='bold'),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=14,face='bold'),
        axis.text.y=element_text(size=14,colour='black'),
        strip.text=element_text(size=12),
        strip.background = element_rect(fill='gray84',colour='black'),
        legend.title=element_text(size=14,face='bold'),
        legend.position = 'right',
        legend.key=element_rect(colour='black',fill='white'),
        legend.spacing.x = unit(0.3, 'cm'),
        panel.background=element_rect(fill='white'),
        panel.border=element_blank(),
        axis.line = element_line(colour = "black"))
        
# Perform statistical analysis
for (i in levels(data$Var2)){
  df<-data[data$Var2==i,]
  print(paste0('Test result for ',i))
  
  var<-var.test(value~disease, df)
  var<-var$p.value
  resi<-shapiro.test(df[df$disease=='Emphysema',]$value)
  resi1<-resi$p.value
  resi<-shapiro.test(df[df$disease=='Fibrosis',]$value)
  resi2<-resi$p.value
  
  if(resi1 | resi2 <= 0.05){
    print(wilcox.test(value ~ disease, data=df))
  }else if(var <= 0.05){
    print(t.test(value ~ disease, data=df, var.equal=FALSE))
  } else{
    print(t.test(value ~ disease, data=df, var.equal=TRUE))
  }
}
