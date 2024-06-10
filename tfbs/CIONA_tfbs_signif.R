
library(MASS)
library(gdata)

setwd('/Users/copley/Dropbox/jobs/emergence/new/tfbs')

tf_p<-function(tf_list,tf_matrix){
   p_vals = c()
   n_tfs = length(tf_list)
   for (i in 1:n_tfs){
      yyy<-tf_list[i]
      fit<-glm(as.formula(paste(yyy,'~ 1 + real')),family='poisson',data=tf_matrix)
#      fit<-glm.nb(as.formula(paste(yyy,'~ 1 + real')),data=tf_matrix)
      p_val<-anova(fit,test="Chisq")$P[2]
#      print(paste(yyy,p_val))
      if(fit$coefficients[2] < 0){ p_val<- 1.0 }
      p_vals<-c(p_vals,p_val)
   }
   p_adj<-p.adjust(p_vals,method='fdr')
   res<-data.frame(p_vals,p_adj)
   row.names(res)<-tf_list
   res
}

#t<-read.table('combined_ciona_tf_matrix_1000up_dust_level_15.txt',header=T)
t<-read.table('combined_ciona_tf_matrix_1000up_no_rc_nseg.txt',header=T)
#t<-read.table('blah.txt',header=T)

# https://stackoverflow.com/questions/23845644/floor-values-in-a-dataframe-in-r
#max = 3
#t[t > max] <- 0

tf_list<-colnames(t)

#tf_list<-c('NFKB2_0')
#tf_list<-c('OTX1_1')

comparisons<-c('aN_vs_aE_T.txt','aE_vs_aN_T.txt',
							 'AN_vs_AM_T.txt','AM_vs_AN_T.txt')

#comparisons<-c('AN_vs_AM_T.txt')

files = file.path('..','results',comparisons)

for(deseq_res in files){

   print(paste('Doing:',deseq_res)) 
	up_list<-deseq_res
   out_file = gsub('_T.txt','_all_tfbs.txt',up_list)

   y<-read.table(up_list,header=T)
   gene_names = strsplit(rownames(y),'_')
   targets<-unlist(lapply(gene_names,'[[',1))
   t['real']<-as.numeric(rownames(t) %in% targets)

   print("Doing real analysis")
   res<-tf_p(tf_list,t)

   signif_tfs<-res[res$p_adj < 0.05,]
   signif_tfs_lst<-rownames(signif_tfs)

   for (j in 1:100){
      print(paste0('control:',j))
      t['real']<-sample(t[,'real'])
      cntl_res<-tf_p(signif_tfs_lst,t)
      cntl_name<-paste0('cntl.',j)
      signif_tfs[cntl_name]<-cntl_res[,'p_vals']
   }

   min_cntrl_p<-apply(signif_tfs[,3:ncol(signif_tfs)], 1, FUN=min)
   final_res<-signif_tfs[1:2]
   final_res['filter']<-as.vector(signif_tfs['p_vals'] < min_cntrl_p)

   resOrdered<-final_res[order(final_res$p_adj),]
   real_matrix<-t[rownames(t) %in% targets,]
   resOrdered['real_counts']<-colSums(real_matrix)[row.names(resOrdered)]
   resOrdered<-resOrdered[c('real_counts','p_vals','p_adj','filter')]
   
   write.fwf(resOrdered,out_file,rownames=T)
#   write.table(resOrdered,out_file,quote=F)
}
