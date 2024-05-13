
library(Seurat)
library(ggplot2)
library(ape)

setwd('./')

raw_counts<-read.table('../counts/raw_counts_w_gene_names.txt',header=T)

batch_design<-read.table('../analysis/batch_design.txt',header=T)

csc<-CreateSeuratObject(counts = raw_counts,min.cells=3,min.features = 200,names.field=1)

csc <- NormalizeData(csc, normalization.method = "LogNormalize", scale.factor = 1000000)

csc <- FindVariableFeatures(csc, selection.method = "vst", nfeatures = 1000)
all.genes <- rownames(csc)
csc <- ScaleData(csc, features = all.genes)
csc <- RunPCA(csc, features = VariableFeatures(object = csc))

Idents(csc)<-batch_design[rownames(umap_dat),]$condition
csc<-BuildClusterTree(csc,dims=1:12,reduction='pca')
PlotClusterTree(csc,show.node.label=False,direction = 'rightwards')
