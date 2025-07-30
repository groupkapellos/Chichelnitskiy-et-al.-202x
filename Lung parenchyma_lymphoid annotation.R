# Remove unneeded objects 
rm(list=ls())

# Load packages
library(readxl)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# Load signatures
features=c('TMSB4X','CD7','IL32','RAD51D','CCL5','ZNF683','KLRD1','PFN1','CCR5','ACTB','NDUFA5','CD27','CD52','NKG7','GZMA','ITGA1','CXCR6','RUNX3',                          # TRM
           'MLLT3','CACNA1A','PDZD3','PIGT','ARSK','HSPA4','SLC25A37','COL6A3','SUV39H2','TTPAL','FAM107A','CDH6','YIPF4','LCMT2','BTBD9','SOX11','RIC3','DOCK7','TNFSF13B',  # TRM-like
           'BTG1','PABPC1','RPL13','CD47','RPL13A','IL7R','JAK3','CD99')                                                                                                      # EM

Snyder 2019 Sci Immunol

features=c('KLRB1','IL7R','SBF2','FOS','LTB','NCR3','CXCR4','CD69','BHLHE40','CD40LG','DPP4','IRF1','GNLY','PRF1',      # MAIT
           'TRDC','TRGC1',                                                                                              # gdT
           'LEF1','ZNF683','XCL1','SELL','TCF7','CXCR3','FCRL3','IL6R','NOSIP','CCR7','FOXP1',                          # memory
           'GZMM','GZMK','IFNG','HLA-DRA','HLA-DRB1','CCL5',                                                            # TEM
           'GZMH','FGFBP2','GZMB','ZEB2','HLA-DRB5','CX3CR1','GZMA','BATF')                                             # memory

features=c('TBX21','NKG7','GZMA','EOMES','SLAMF7','GZMK','PRF1','FCRL3','MAF','CCR5','GZMH','TIGIT','TOX','PDCD1','FGFBP2','CX3CR1','IFNG','ZEB2') # exhausted

features=c('LEF1','SELL','CCR7','BACH2','SATB1','NOSIP','IL6R','RUNX2','GZMB','RORA','XBP1','RORC','SOCS1','KLRD1','DPP4','GNLY','IL26','IL23R','IFIT2','IFIT3','IRF8','IL2','ID3','GATA3',                                                                 # stem
           'TOX','TOX2','TIGIT','PDCD1','NFATC2','EOMES','BATF','ZEB2','PRF1','GZMA','TNFSF9','CCL18','CCL4','CCL3','CD84','CD74','SLAMF6','MAF','CD80','IKZF3','CD44','FAS','FCRL3','KLF3','SMAD3','GZMH','CXCR3','GZMK','SOX5','CXCL13','CCR5','TNFRSF9') # exhausted

Galletti 2020 Nat Immunol

features=c('ID2','RUNX3','PRDM1','TOX',                       # terminal differentiation
           'GZMA','GZMB','GZMH','PRF1','GNLY','NKG7',         # cytotoxic
           'ICOS','TNFRSF4','TNFRSF18',                       # co-stimulation
           'CXCL13','CX3CR1','CXCR3',                         # chemokines/chemokine receptors
           'IL2RA','FOXP3','CTLA4','RTKN2','TNFRSF15')        # T reg

features=c('ITGAE','ITGA1','CXCR6','RBPJ','ZNF683','PDCD1','KLRC1','SPRY1',   # TRM
           'S1PR5','KLF3','CX3CR1','S1PR1','KLRG1','EOMES',                   # TRM-like
           'TCF7', 'CCR7','SELL')                                             # TCM

Clarke 2019 JEM

features=c('CD3D','CD3E','CD3G','CD4','CD8A','CD8B',
           'CD40LG','CXCR6',                                                        # CD4 EM
           'ITGA1','ITGAE',                                                         # CD4 TRM
           'LTB','LEF1','CD28','KLF2','SELL',                                       # TCM
           'GZMH','KLRG1',                                                          # TEMRA
           'GZMK','EOMES','DTHD1','CRTAM','LYST','TNIP3','PDCD1',                   # CD8 TRM/EM
           'AREG','SOX4','KIT','TNFRSF18',                                          # ILC
           'DUSP2','TRBV6-2','SLC4A10','IFNGR1','IL7R','KLRB1','CCR6',              # MAIT cells
           'NKG7','FGFBP2','KIR2DL1','KLRF1','GNLY','FCER1G',                       # CD16hi NK cells
           'S1PR1','IL32',                                                          # NKT cells
           'KLRC3','ITGAX','KIR3DX1','CD247','ITGAD','PFN1','GZMA','GZMB',          # CD11d- NK cells
           'IL2RB','NCR1','CCL3','CD7','CMC1','XCL1','CCL4','CCL5',                 # CD56bri NK cells
           'KLRC2','TRDC','TRG-AS1','TRGC1','KIR2DL4',                              # gd
           'FOXP3','CTLA4','CCR4','IL2RA','TNFRSF4','TIGIT','ZNF683')               # T reg

Madissoon 2022 Nat Genet

features=c('TCF7','ID3','CCR7','AQP3','SELL',                                                            # naive
           'IL4R','STAT1','MAL','SOCS1','IL2','ODC1','PSAT1','WARS','PYCR1','TNF','MIR155HG','NME1',     # activated
           'FOXP3','RGS1','LGALS3','IL2RA','CTLA4','TIGIT','TNFRSF4','TNFRSF18',                         # Tregs
           'LGALS1','GZMK','CXCR6','ANXA2','ITGA1',                                                      # CD8 TRM/EM
           'GZMB','IFNG','XCL1','CCL4','CCL3',                                                           # CD8 TRM/EM activated
           'GZMH','NKG7','MYO1F','CCL5','GNLY',                                                          # CD8 TEMRA
           'IFIT3','KLRD1','XCL2','HOPX','RCAN2','PRF1')                                                 # CD8 TEMRA activated

Szabo 2019 Nat Commun

features=c('KLRF1','CCL3','KLRB1','FCGR3A','SPON2','FGFBP2','PLAC8','GNLY','PRF1','GZMB','GZMH',                                                                                                                        # NK                                                        
           'FGFBP2','PLAC8','GNLY','PRF1','GZMB','TRAC','CD3D','GZMH','IL7R',                                                                                                                                           # T NK-like
           'CCL3','KLRB1','GNLY','PRF1','GZMB','TRAC','CD3D','CD8A','CD8B','GZMH','GZMK','IL7R',                                                                                                                        # T GZMK
           'CCL3','KLRB1','PRF1','GZMB','TRAC','CD3D','CD8A','CD8B','GZMH','GZMK','IL7R','HLA-DRA','HLA-DRB1',                                                                                                          # CD8 TRM
           'CCL3','KLRB1','PRF1','GZMB','TRAC','CD3D','CD8A','CD8B','GZMH','CD27','IFNG','ITGA1','ITGAE','PDCD1','TOX','HAVCR2','LAG3','CXCL13','ENTPD1','TIGIT','TNFRSF9','TNFRSF18','CTLA4','HLA-DRA','HLA-DRB1',     # T activated
           'KLRB1','GZMB','TRAC','CD3D','IL7R','HLA-DRA','HLA-DRB1',                                                                                                                                                    # CD4 TRM
           'KLRB1','TRAC','CD3D','IL7R',                                                                                                                                                                                # naive/TCM
           'TRAC','CD3D','IL7R','CD27','TIGIT','TNFRSF9','TNFRSF18','CTLA4','ICOS','MAGEH1','IL2RA','FOXP3','TNFRSF4','HLA-DRA','HLA-DRB1',                                                                             # Tregs
           'TRAC','CD3D','HLA-DRA','HLA-DRB1','MKI67','STMN1')                                                                                                                                                          # Cycling

features=c('MS4A1','HLA-DRB1','CD69','LTB','GPR183','CD83','BTG1','CXCR4','ACTB',          # B                                                         
           'MZB1','XBP1','SSR4','IGHD','JCHAIN',                                           # IgD plasma
           'IGHM',                                                                         # IgM plasma
           'IGHA1','IGHA2',                                                                # IgA plasma
           'IGHG1','IGHG2','IGHG3','IGHG4')                                                # IgG plasma

Leader 2021 Cancer Cell

features=c('IFI6','ISG15','IFI44L','XAF1','MX1',
           'MKI67','LGALS1','STMN1','TUBA1B','HMGB2',
           'CCL4','CCL4L2','CCL3','CCL3L3','IFNG',
           'HLA-DPA1','HLA-DRA','HLA-DRB1','HLA-DPB1','CD74',
           'PRF1','GZMB','FCGR3A',
           'NCAM1','SELL','GZMK','TCF7','XCL1')     

DotPlot(nks, features=features, dot.min=0.1, cols='RdBu') + RotatedAxis()

# Gene signature heatmaps
test<-nks@assays$RNA@data  
test<-test[rowSums(as.matrix(test))>0,]

tmp2<-data.frame(rownames(test))
colnames(tmp2)<-'genes'

for (i in levels(nks$integrated_snn_res.0.3)) {
  tmp<-test[,colnames(test) %in% rownames(nks@meta.data[nks$integrated_snn_res.0.3==i,])]
  tmp<-as.data.frame(tmp)
  tmp<-as.data.frame(rowMeans(tmp))
  colnames(tmp)<-i
  tmp2<-cbind(tmp2,tmp)
}  
tmp2$genes<-NULL

features=c('TMSB4X','IL32','CCL5','CD52','ITGA1','RBPJ','RGS1','IFNGR1','FOS','CXCR4','ID2',
           'CD69','BTG1','PABPC1','RPL13','RPL13A','IL7R',
           'SELL','CCR7','IL2RA','FOXP1',
           'LTB',
           'STMN1','MKI67','LGALS1',
           'PFN1','KLRD1','PRF1','NKG7','GZMH','RUNX3','FGFBP2','LYST','ZEB2','CD47','NDUFA5','IKZF3',
           'CD99','KLRB1','KLRF1','PLAC8','SLAMF7','KLF3','S1PR5','MYO1F','CRTAM','CD7','DUSP2',
           'DTHD1','RORA','CCL3','GZMB','CCL4','GNLY','GZMA','CD44','CD74','HLA-DRA','HLA-DRB1','HOPX','LGALS3','STAT1','IFNG','BTBD9','BHLHE40','XBP1','ANXA2')

# Snyder, Galletti, Clarke, Madissoon, Szabo, Leader

features=c('MS4A1','HLA-DRB1','CD69','LTB','GPR183','CD83','BTG1','CXCR4','ACTB',          # B                         
           'IGHM',                                                                         # IgM plasma
           'IGHG1','IGHG2','IGHG3','IGHG4',                                                # IgG plasma
           'MZB1','XBP1','SSR4','JCHAIN',                                                  # plasma
           'IGHA1','IGHA2')                                                                # IgA plasma
                                                                                                     
test2<-tmp2[rownames(tmp2) %in% unique(features),]
test2<-test2[unique(features),]
test2<-test2[!grepl('NA', rownames(test2)),]

pheatmap(test2,
         clustering_distance_rows='correlation',
         clustering_distance_cols='correlation',
         scale='row',
         show_colnames=T,
         cluster_rows=F,
         cluster_cols=F,
         cellheight=8,
         cellwidth=12,
         breaks=seq(-2, 2, by=0.05),
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(length(seq(-2, 2, by=0.05))))