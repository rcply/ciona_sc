
library('DESeq2')
library('ashr')

setwd('./')

#
# read in full count table
#

filtering = 'T'
p_adj_cut<-0.01
lfc_cut<-0.0

count_table<-read.table('counts/raw_counts_w_gene_names.txt',header=T)
outlier_cells<-c('A61l5','A87r5')
remove_list<-which(colnames(count_table) %in% outlier_cells)
if(length(remove_list) > 0)
   count_table<-count_table[,-remove_list]

comparisons<-read.table('comparisons_table.txt',header=T,sep='\t')
#comparisons<-read.table('test_table.txt',header=T,sep='\t')

for (i in 1:nrow(comparisons)){

   cell_a <- comparisons[i,'cell_a']
   cell_b <- comparisons[i,'cell_b']
   all_cell_a <- grep(cell_a,colnames(count_table),value=T)
   all_cell_b <- grep(cell_b,colnames(count_table),value=T)
   all_cells <- c(all_cell_a,all_cell_b)
   print(all_cells)
   cell_set_counts<-count_table[,all_cells]
   batch_design<-read.table(comparisons[i,'design'],header=T)
   cell_set_design<-batch_design[drop=F,rownames(batch_design) %in% 
                           colnames(cell_set_counts),]
   dds<-DESeqDataSetFromMatrix(countData = cell_set_counts,
                            colData = cell_set_design[colnames(cell_set_counts),], 
                            design = as.formula(comparisons[i,'formula']))
   dds<-DESeq(dds)
   
   res<-results(dds,alpha=p_adj_cut,lfcThreshold=lfc_cut,
   						 contrast=c('condition',cell_a,cell_b),
   						 independentFiltering = filtering
   						 )
   
   resLFC <- lfcShrink(dds,contrast=c('condition',cell_a,cell_b),res=res,type='ashr')
   #resLFC <- lfcShrink(dds,coef=2,type='apeglm')

   resOrdered<-resLFC[order(resLFC$padj),]
   y<-na.omit(resOrdered)
   yy<-y[y$padj < p_adj_cut,]
   
   out_fn = paste0(cell_a,'_vs_',cell_b,'_',filtering,'.txt')
   a_vs_b <-yy[yy$log2FoldChange > 0,]
   write.table(a_vs_b,out_fn,quote=F)
   out_fn = paste0(cell_b,'_vs_',cell_a,'_',filtering,'.txt')
   b_vs_a <-yy[yy$log2FoldChange < 0,]
   write.table(b_vs_a,out_fn,quote=F)
   
}
