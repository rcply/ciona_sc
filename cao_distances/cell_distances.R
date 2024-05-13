library(pheatmap)
library(RColorBrewer)

setwd('./')

t<-read.table('cells_no_Pb.txt',header=F,row.names=1)
key_table<-read.table('key_to_name_no_Pb.txt',sep='\t',header=T)

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# dice
library(philentropy)
dist_obj<-distance(t,'dice',use.row.names = T)
dist_mat<-as.matrix(dist_obj)
cell_names = arrange(key_table,sapply(key,function(y) which(y==row.names(dist_mat))))$name
cell_lineage = substr(row.names(dist_mat),1,1)
pheatmap(dist_mat,col=colors,treeheight_col=0,fontsize=8,
         labels_col = cell_lineage,labels_row = cell_names,angle_col=0)
