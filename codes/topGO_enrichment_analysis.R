library(topGO)
library(tidyr)
library(dplyr)
library(tibble)

readdir_cluster="E:/Grad study/lab/data/soft kmeans/gene selection/ygob_libnorm_lbt005_cluster_result/cluster file/k="
filelist=list.files(path="E:/Grad study/lab/data/soft kmeans/gene selection/ygob_libnorm_lbt005_cluster_result/cluster file/k=3")
savedir="E:/Grad study/lab/data/soft kmeans/gene selection/ygob_libnorm_lbt005_cluster_result/topgo_result/raw pre 241031"
gene_list=read.csv("E:/Grad study/lab/data/soft kmeans/gene selection/ygob_lbt001_all_4_species_with_orth.csv")[,1:3]
k_list=c(3,4,5,7,9)
# this allows you to select the needed w to be plotted
reward=c(0,0.05,0.1,0.2,0.5,0.6,0.7,0.8,0.9,1,1.5) # change yourself
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
    assign(paste0(k,"_",a),read.csv(paste0(readdir_cluster,k,"/",a)))
    
  }
}


for(k in k_list){
  for(a in filelist){
    assign(paste0(k,"_",a),merge(get(paste0(k,"_",a)),gene_list,by=c("orth_group","num")))
    assign(paste0(k,"_",a),get(paste0(k,"_",a))[,c(4,1:3)])
    assign(paste0(k,"_",a),get(paste0(k,"_",a))[order(get(paste0(k,"_",a))$group),])
  }
}


## take k=3, w=0.1 as example
###################
k=4
a="w=0.csv"
readdir_GO="E:/Grad study/lab/data/GO terms/annotated GO/gene_GO_"
gene2GO_sc=read.csv(paste0(readdir_GO,"sc",".csv"))
gene2GO_ca=read.csv(paste0(readdir_GO,"ca",".csv"))
gene2GO_cg=read.csv(paste0(readdir_GO,"cg",".csv"))
gene2GO_kl=read.csv(paste0(readdir_GO,"kl",".csv"))


gene2GO_df=rbind(gene2GO_ca,gene2GO_cg,gene2GO_sc)

gene2GO=gene2GO_df[,c("Gene","GO")] %>%
  group_by(Gene) %>%
  summarize(GO = list(GO)) %>%
  deframe()

allGenes <- names(gene2GO)

# Create a named vector with 1 for the genes of interest and 0 for the rest
# First analyze cluster 1
for(g in 1:k){
  geneList_g1 <- factor(as.integer(allGenes %in% get(paste0(k,"_",a))$Geneid[get(paste0(k,"_",a))$group==g]))  # Adjust for your actual gene list
  names(geneList_g1) <- allGenes
  for(go in c("CC","BP","MF")){
    GOdata_g1_BP <- new("topGOdata",
                        ontology = go,  # Choose "BP" for Biological Process, or "MF" or "CC"
                        allGenes = geneList_g1,
                        annot = annFUN.gene2GO,
                        gene2GO = gene2GO)
    resultFisher <- runTest(GOdata_g1_BP, algorithm = "classic", statistic = "fisher")
    # Define GOdata object for analysis
    GOtable_g1_BP=GenTable(GOdata_g1_BP, classicFisher = resultFisher,topNodes=10)
    # View the top results
    #GenTable(GOdata_g1_BP, classicFisher = resultFisher,topNodes=10)
    write.csv(GOtable_g1_BP,paste0(savedir,"/k=",k,"_",strsplit(a,".csv")[[1]],"_g",g,"_",go,".csv"))
  }
}






############################

## run topGO without Kl
############################

readdir_GO="E:/Grad study/lab/data/GO terms/annotated GO/gene_GO_"
gene2GO_sc=read.csv(paste0(readdir_GO,"sc",".csv"))
gene2GO_ca=read.csv(paste0(readdir_GO,"ca",".csv"))
gene2GO_cg=read.csv(paste0(readdir_GO,"cg",".csv"))
gene2GO_kl=read.csv(paste0(readdir_GO,"kl",".csv"))


gene2GO_df=rbind(gene2GO_ca,gene2GO_cg,gene2GO_sc)

gene2GO=gene2GO_df[,c("Gene","GO")] %>%
  group_by(Gene) %>%
  summarize(GO = list(GO)) %>%
  deframe()

allGenes <- names(gene2GO)
for(k in 3:12){
  for(a in filelist){
    w=strsplit(a,".csv")[[1]][1]
    for(g in 1:k){
      for(go in c("BP","MF","CC")){
        geneList <- factor(as.integer(allGenes %in% get(paste0(k,"_",a))$Geneid[get(paste0(k,"_",a))$group==g]))  # Adjust for your actual gene list
        names(geneList) <- allGenes
        
        # Define GOdata object for analysis
        GOdata <- new("topGOdata",
                            ontology = go,  # Choose "BP" for Biological Process, or "MF" or "CC"
                            allGenes = geneList,
                            annot = annFUN.gene2GO,
                            gene2GO = gene2GO)
        resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
        

        GOtable=GenTable(GOdata, classicFisher = resultFisher,topNodes=300)
        # 
        # savedir_1=paste0(savedir,"/topGONodes/k=",k)
        # savedir_2=paste0(savedir_1,"/",w)
        # if(!dir.exists(savedir_1)){
        #   dir.create(savedir_1)
        # }
        # if(!dir.exists(savedir_2)){
        #   dir.create(savedir_2)
        # }
        savedir_3=paste0(savedir,"/original_Fisher_result/k=",k)
        savedir_4=paste0(savedir_3,"/",w)
        if(!dir.exists(savedir_3)){
          dir.create(savedir_3)
        }
        if(!dir.exists(savedir_4)){
          dir.create(savedir_4)
        }
        savedir_5=paste0(savedir,"/topGO_data/k=",k)
        savedir_6=paste0(savedir_5,"/",w)
        if(!dir.exists(savedir_5)){
          dir.create(savedir_5)
        }
        if(!dir.exists(savedir_6)){
          dir.create(savedir_6)
        }
        #write.csv(GOtable,paste0(savedir_2,"/g",g,"_",go,".csv"),row.names = FALSE)
        save(resultFisher, file=paste0(savedir_4,"/g",g,"_",go,'.rda'))
        save(GOdata, file=paste0(savedir_6,"/g",g,"_",go,'.rda'))
      }
    }
  }
}
