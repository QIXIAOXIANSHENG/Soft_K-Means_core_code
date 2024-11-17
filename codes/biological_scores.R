library(topGO)

# tryout codes, don't run
###########################################
fisher_path <- "D:/Grad study/lab/data/soft kmeans/gene selection/ygob_lbt001_cluster_result/topgo_result/without Kl/original_Fisher_result/k=3"
godata_path <- "D:/Grad study/lab/data/soft kmeans/gene selection/ygob_lbt001_cluster_result/topgo_result/without Kl/topGO_data/k=3"
# leave the rest untouched unless you have different folder names.
items <- list.files(fisher_path)
for(w in items){
  for(g in 1:3){
    load(paste0(fisher_path,"/",w,"/g",g,"_BP.rda"))
    load(paste0(godata_path,"/",w,"/g",g,"_BP.rda"))
    df=GenTable(GOdata,classicFisher = resultFisher,topNodes=length(GOdata@graph@nodes))[,c(1,3,4)]
    df$score=sort(score(resultFisher))
    assign(paste0(w,"g",g),df)
    
  }
}
for(w in items){
  for(g in 1:3){
    df=get(paste0(w,"g",g))
    assign(paste0(w,"g",g),df[df$Annotated>30 & df$score<0.05,])
  }
}


for(w in items){
  assign(paste0("T_acat_uniweight",w),0)
  for(g in 1:3){
    
    assign(paste0("T_acat_uniweight",w),
           get(paste0("T_acat_uniweight",w))+
             mean(tan((0.5-get(paste0(w,"g",g))$score)*pi)))
  }
}

for(w in items){
  assign(paste0("p_acat_uniweight",w),0.5-atan(get(paste0("T_acat_uniweight",w)))/pi)
  cat(w,":",get(paste0("p_acat_uniweight",w)),"\n")
}



for(w in items){
  assign(paste0("T_acat_weight_gaussian_self",w),0)
  weight=c()
  for(g in 1:3){
    bw=(max(log(get(paste0(w,"g",g))$score))-min(log(get(paste0(w,"g",g))$score)))/nrow(get(paste0(w,"g",g)))
    kde_result <- density(log(get(paste0(w,"g",g))$score), kernel = "gaussian", bw = bw) 
    weight_g=nrow(get(paste0(w,"g",g)))*approx(kde_result$x, kde_result$y, xout = log(get(paste0(w,"g",g))$score))$y
    weight= append(weight,weight_g)
    t_g=sum(weight_g*tan((0.5-get(paste0(w,"g",g))$score)*pi))
    assign(paste0("T_acat_weight_gaussian_self",w),
           get(paste0("T_acat_weight_gaussian_self",w))+
             t_g)
  }
  p=0.5-atan(get(paste0("T_acat_weight_gaussian_self",w))/sum(weight))/pi
  assign(paste0("p_acat_weight_gaussian_self",w),p)
}

bw_choose=c(0.2,0.5,1,1.5,2,3)
for(bw in bw_choose){
  for(w in items){
    assign(paste0("T_acat_weight_gaussian_",bw,w),0)
    weight=c()
    for(g in 1:3){
      kde_result <- density(log(get(paste0(w,"g",g))$score), kernel = "gaussian", bw = bw) 
      weight_g=nrow(get(paste0(w,"g",g)))*approx(kde_result$x, kde_result$y, xout = log(get(paste0(w,"g",g))$score))$y
      weight= append(weight,weight_g)
      t_g=sum(weight_g*tan((0.5-get(paste0(w,"g",g))$score)*pi))
      assign(paste0("T_acat_weight_gaussian_",bw,w),
             get(paste0("T_acat_weight_gaussian_",bw,w))+
               t_g)
    }
    p=0.5-atan(get(paste0("T_acat_weight_gaussian_",bw,w))/sum(weight))/pi
    assign(paste0("p_acat_weight_gaussian_",bw,w),p)
  }
}


for(w in items){
  cat(w,":",get(paste0("p_acat_weight_gaussian_self",w)),"\n")
}

df_k=as.data.frame(matrix(NA,1,length(bw_choose)+2))
colnames(df_k)=c("w",paste0("score_",bw_choose),"score_self")
for(w in items){
  
  go=c()
  for(g in 1:3){
    go=append(go,get(paste0(w,"g",g))$GO.ID)
  }
  n=length(unique(go))
  if(w=="w=0"){
    n0=n
  }
  b=n0/n
  vec=c()
  for(bw in c(bw_choose,"self")){
  a=get(paste0("p_acat_weight_gaussian_",bw,w))/get(paste0("p_acat_weight_gaussian_",bw,"w=0"))
  vec=append(vec,a+b)
  }
  df_temp=as.data.frame(matrix(c(strsplit(w,"=")[[1]][2],vec),1))
  colnames(df_temp)=c("w",paste0("score_",bw_choose),"score_self")
  
  #smaller the better
  if(w=="w=0"){
    df_k=df_temp
  }
  else{
    df_k=rbind(df_k,df_temp)
  }
}


###########################################


# true codes
###########################################
fisher_path <- "D:/Grad study/lab/data/soft kmeans/gene selection/ygob_lbt001_cluster_result/topgo_result/without Kl/original_Fisher_result/k="
godata_path <- "D:/Grad study/lab/data/soft kmeans/gene selection/ygob_lbt001_cluster_result/topgo_result/without Kl/topGO_data/k="
# leave the rest untouched unless you have different folder names.
bw_choose=c(0.2,0.5,1,1.5,2,3,"self")
k=3
items <- list.files(paste0(fisher_path,k))
df_k=data.frame(matrix(NA,(12-3+1)*length(items),length(bw_choose)+2))
colnames(df_k)=c("k","w",paste0("score_",bw_choose))
df_k$k=rep(c(3:12),each=length(items))
df_k$w=rep(unlist(lapply(items,function(x){strsplit(x,"=")[[1]][2]})),(12-3+1))
for(k in 3:12){
  for(bw in bw_choose){
    score=c()
    for(w in items){
      weight=c()
      T_acat=0
      go=c()
      for(g in 1:k){
        load(paste0(fisher_path,k,"/",w,"/g",g,"_BP.rda"))
        load(paste0(godata_path,k,"/",w,"/g",g,"_BP.rda"))
        df=GenTable(GOdata,classicFisher = resultFisher,topNodes=length(GOdata@graph@nodes))[,c(1,3,4)]
        df$score=sort(score(resultFisher))
        wg=df[df$Annotated>30 & df$score<0.05,]
        go=append(go,wg$GO.ID)

        if(bw!="self"){
          kde_result <- density(log(wg$score), kernel = "gaussian", bw = as.numeric(bw))
        }
        else{
          bw_1=(max(log(wg$score))-min(log(wg$score)))/nrow(wg)
          kde_result <- density(log(wg$score), kernel = "gaussian", bw = bw_1)
        }
        weight_g=nrow(wg)*approx(kde_result$x, kde_result$y, xout = log(wg$score))$y
        weight= append(weight,weight_g)
        t_g=sum(weight_g*tan((0.5-wg$score)*pi))
        T_acat=T_acat+t_g
      }
      p=0.5-atan(T_acat/sum(weight))/pi
      if(w=="w=0"){
        p0=p
      }
      n=length(unique(go))
      if(w=="w=0"){
        n0=n
      }
      b=n0/n
      score=append(score,p/p0+b)
    }
    df_k[df_k$k==k,paste0("score_",bw)]=score
  }
}

write.csv(df_k,"D:/Grad study/lab/data/soft kmeans/gene selection/ygob_lbt001_cluster_result/topgo_result/without Kl/clustering scores/gaussian_bio_score.csv",index=FALSE)
