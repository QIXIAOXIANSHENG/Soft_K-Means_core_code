######
################

orth=read.table("D:/Grad study/lab/data/soft kmeans/ortho information from ygob/20240206-4sps-ygob-cgob-one-gene-per-sps.tsv",sep="\t",header=TRUE)
data=read.csv("D:/Grad study/lab/data/soft kmeans/gene selection/ygob_libnorm_lbt005_all_4_species.csv")

###### get all unique genes
a=unique(data$Geneid)
del=c()
for(i in a){
  b=which(data$Geneid==i)
  if(length(b)>1){
    del=append(del,b[-1])
  }
}
data=data[-del,]

###### get orthogroup label
data$orth=rep(NA,nrow(data))
data$num=rep(NA,nrow(data))


for(i in 1:nrow(orth)){
  list=as.vector(orth[i,-7][which(!is.na(orth[i,-7]))])
  for (j in 1:length(list)){
    
    data$orth[which(list[[j]]==data$Geneid)]=i
    data$num[which(list[[j]]==data$Geneid)]=j
  }
}
data_omit=na.omit(data)
data_omit_1=data_omit[,c(1,11,12,2:10)]


for(i in unique(data_omit_1$orth)){
  if(length(data_omit_1[data_omit_1$orth==i,"num"])>6){
    print(i)
  }
}
#data_omit_1=data_omit_1[-which(data_omit_1$orth==2968),]
#data_omit_1=data_omit_1[-which(data_omit_1$orth==3261),]
#data_omit_1=data_omit_1[-which(data_omit_1$orth==454),]
names(data_omit_1)[2]="orth_group"
del=c()
for(i in unique(data_omit_1$orth_group)){
  if(length(data_omit_1[data_omit_1$orth_group==i,"num"])!=max(data_omit_1[data_omit_1$orth_group==i,"num"])){
   del=append(del,i)
  }
}
for(i in del){
  data_omit_1[data_omit_1$orth_group==i,"num"]=1:length(data_omit_1[data_omit_1$orth_group==i,"num"])
}

a=apply(data_omit_1,1,function(x){
  paste0(x[2],"_",x[3])
})
length(unique(data_omit_1$orth_group))

data_omit_1=data_omit_1[order(data_omit_1$orth_group,data_omit_1$num),]

write.csv(data_omit_1,
          "D:/Grad study/lab/data/soft kmeans/gene selection/library_norm/ygob_libnorm_lbt005_all_4_species_with_orth.csv",
          row.names = FALSE)




