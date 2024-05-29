
library('DESeq2')
library('ggplot2')

setwd('./')

count_table<-read.table('counts/raw_counts_w_gene_names.txt',header=T)

batch_design<-read.table('batch_design_for_figures.txt',header=T)

count_table<-count_table[,rownames(batch_design)]

outlier_cells<-c('A61l5','A87r5')
remove_list<-which(colnames(count_table) %in% outlier_cells)
if(length(remove_list) > 0)
   count_table<-count_table[,-remove_list]         

cell_set_counts<-count_table

cell_set_design<-batch_design[drop=F,rownames(batch_design) %in% 
                           colnames(cell_set_counts),]

dds<-DESeqDataSetFromMatrix(countData = cell_set_counts, colData = cell_set_design, design = ~ condition) 

dds<-DESeq(dds)

vsd<-vst(dds,blind=T)

# basics of plot with all data points
pcaData<-plotPCA(vsd,returnData=T)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape = condition)) + 
    geom_point(size =2.5,stroke=1) + 
   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
   ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
     scale_shape_manual(values=seq(0,18)) +
     coord_fixed() + theme(title=element_text(size=18), legend.text=element_text(size=18))
