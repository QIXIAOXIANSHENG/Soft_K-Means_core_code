ca_raw=read.csv("D:\\Grad study\\lab\\data\\soft kmeans\\Ca\\Ca_rawCounts.csv")
cg_raw=read.csv("D:\\Grad study\\lab\\data\\soft kmeans\\Cg\\Cg_rawCounts.csv")
sc_raw=read.csv("D:\\Grad study\\lab\\data\\soft kmeans\\Sc\\Sc_rawCounts.csv")
kl_raw=read.csv("D:\\Grad study\\lab\\data\\soft kmeans\\Kl\\Kl_rawCounts.csv")

row.names(ca_raw)=ca_raw$Geneid
row.names(cg_raw)=cg_raw$Geneid
row.names(sc_raw)=sc_raw$Geneid
row.names(kl_raw)=kl_raw$Geneid

ca_raw=ca_raw[,-1]
cg_raw=cg_raw[,-1]
sc_raw=sc_raw[,-1]
kl_raw=kl_raw[,-1]

sc_lib=apply(sc_raw,2,sum)
cg_lib=apply(cg_raw,2,sum)
ca_lib=apply(ca_raw,2,sum)
kl_lib=apply(kl_raw,2,sum)

sc_raw=t(apply(sc_raw,1,function(x){x*1000/sc_lib}))
cg_raw=t(apply(cg_raw,1,function(x){x*1000/cg_lib}))
ca_raw=t(apply(ca_raw,1,function(x){x*1000/ca_lib}))
kl_raw=t(apply(kl_raw,1,function(x){x*1000/kl_lib}))

# DESeq
###############################
# using the first two reads as control
coldata_ca <- data.frame(condition = factor(c(rep('control',2),rep("treat",24)), levels = c('control', 'treat')))
coldata_cg <- data.frame(condition = factor(c(rep('control',2),rep("treat",20)), levels = c('control', 'treat')))
coldata_sc <- data.frame(condition = factor(c(rep('control',2),rep("treat",20)), levels = c('control', 'treat')))
coldata_kl <- data.frame(condition = factor(c(rep('control',2),rep("treat",24)), levels = c('control', 'treat')))

library(DESeq2)
for(n in c("ca","cg","sc","kl")){
 dds=DESeqDataSetFromMatrix(countData = get(paste0(n,"_raw")), colData = get(paste0("coldata_",n)), design= ~condition)
  dds1=DESeq(dds)
  res=results(dds1, contrast = c('condition', 'treat', 'control'))
  res=res[order(res$padj),]
  for(p in c(0.05,0.1,0.15)){
    diff_gene=row.names(subset(res,padj < p))
    resdata=get(paste0(n,"_raw"))[diff_gene,]
    write.csv(resdata,file= paste0("D:\\Grad study\\lab\\data\\soft kmeans\\gene selection\\",n,"_2ndcontrol_padj_",p,".csv"),row.names = TRUE)
  }
}


# dds_ca <- DESeqDataSetFromMatrix(countData = ca_raw, colData = coldata_ca, design= ~condition)
# dds1_ca <- DESeq(dds_ca)
# res_ca <- results(dds1_ca, contrast = c('condition', 'treat', 'control'))
# summary(res_ca)
# # out of 6486 with nonzero total read count
# # adjusted p-value < 0.1
# # LFC > 0 (up)       : 1250, 19%
# # LFC < 0 (down)     : 920, 14%
# # outliers [1]       : 0, 0%
# # low counts [2]     : 0, 0%
# # (mean count < 0)
# # [1] see 'cooksCutoff' argument of ?results
# # [2] see 'independentFiltering' argument of ?results
# res_ca <- res_ca[order(res_ca$padj),]
# # table(res_ca$padj<0.15) 
# # FALSE  TRUE 
# # 4041  2445 
# diff_gene_ca<-row.names(subset(res_ca,padj < 0.05))
# 
# resdata_ca <- ca_raw[diff_gene_ca,]
# write.csv(resdata_ca,file= "D:\\Grad study\\lab\\data\\soft kmeans\\gene selection\\ca_2ndcontrol_padj_005.csv",row.names = TRUE)
###############################

# test for white noise
################################
for( n in c("ca","cg","sc","kl")){
  assign(paste0("lbt_p_",n),c())
  for(i in 1:nrow(get(paste0(n,"_raw")))){
    assign(paste0("lbt_p_",n),append(get(paste0("lbt_p_",n)),
                                     mean(unlist(
                                       lapply(1:10,function(x){
                                         Box.test(ts(as.vector(as.matrix(get(paste0(n,"_raw"))[i,]))), lag = x, type = "Ljung-Box")[3]
                                         })))))
    
  }
}

for( n in c("ca","cg","sc","kl")){
assign(paste0("lbt_",n),get(paste0(n,"_raw"))[which(get(paste0("lbt_p_",n))<0.05),])
}

for( n in c("ca","cg","sc","kl")){
write.csv(get(paste0("lbt_",n)),file= paste0("D:\\Grad study\\lab\\data\\soft kmeans\\gene selection\\results_LBtest\\libnorm_",n,"_lbt_p_005.csv"),row.names = TRUE)
}

################################
