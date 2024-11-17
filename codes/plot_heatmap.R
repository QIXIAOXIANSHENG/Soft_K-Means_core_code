readdir="E:/Grad study/lab/data/soft kmeans/gene selection/ygob_libnorm_lbt005_cluster_result/plot_heatmap_csv file/k="
filelist=list.files(path="E:/Grad study/lab/data/soft kmeans/gene selection/ygob_libnorm_lbt005_cluster_result/plot_heatmap_csv file/k=3")
savedir="E:/Grad study/lab/data/soft kmeans/gene selection/ygob_libnorm_lbt005_cluster_result"
k_list=c(3,4,5,7,9)

# this allows you to select the needed w to be plotted
reward=c(0.6,0.7,0.8,0.9,1,1.5) # change yourself
new_filelist=c()
for(a in filelist){
  w_=strsplit(strsplit(a,".csv")[[1]],"=")[[1]][2]
  if(w_ %in% reward){
    new_filelist=append(new_filelist,a)
  }
}
filelist=new_filelist
####

for(k in k_list){
  for(a in filelist){
    assign(paste0(k,"_",a),read.csv(paste0(readdir,k,"/",a)))
    
  }
}
gene_list=read.csv("E:/Grad study/lab/data/soft kmeans/gene selection/library_norm/ygob_libnorm_lbt005_all_4_species_with_orth.csv")[,1:3]

for(k in k_list){
  for(a in filelist){
    assign(paste0(k,"_",a),merge(get(paste0(k,"_",a)),gene_list,by=c("orth_group","num")))
    assign(paste0(k,"_",a),get(paste0(k,"_",a))[,c(13,1:12)])
    assign(paste0(k,"_",a),get(paste0(k,"_",a))[order(get(paste0(k,"_",a))$group),])
  }
}

library(pheatmap)


for(k in k_list){
  savedir_1=paste0(savedir,"/plot_heatmap file")
  if(!dir.exists(savedir_1)){
    dir.create(savedir_1)
  }
  savedir_k=paste0(savedir_1,"/k=",k)
  if(!dir.exists(savedir_k)){
    dir.create(savedir_k)
  }
  for(a in filelist){
    data_matrix=as.matrix(get(paste0(k,"_",a))[,-c(1,2,3,4)])
    rownames(data_matrix) <- get(paste0(k,"_",a))$Geneid
    w=strsplit(a,split=".csv")[[1]][1]
    anno <- data.frame(group = factor(get(paste0(k, "_", a))$group), 
                       row.names = rownames(data_matrix))
    
    # Create a color palette for the unique groups
    newCols <- colorRampPalette(grDevices::rainbow(length(unique(anno$group))))
    annoCol <- newCols(length(unique(anno$group)))
    
    # Assign names to the color palette based on the unique levels in the factor
    names(annoCol) <- levels(anno$group)
    
    # Create the annotation_colors list with the correct structure
    annoCol <- list(group = annoCol)
    
    pheatmap(data_matrix,
             annotation_row = anno,
             annotation_colors = annoCol,
             cluster_row = FALSE, cluster_cols = FALSE,
             filename=paste0(savedir_k,"/heatmap ",w,".png"),
             height=10,width = 5.51,
             main=paste0("heatmap ",w),show_rownames=FALSE)
    
  }
}
