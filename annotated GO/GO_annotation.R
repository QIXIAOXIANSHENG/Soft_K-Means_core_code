library(gprofiler2)
library(tidyr)
library(dplyr)
library(tibble)

gene_list=read.csv("D:/Grad study/lab/data/soft kmeans/gene selection/ygob_lbt001_all_4_species_with_orth.csv")[,1]
kl_genelist=grep("^KLL", gene_list, value = TRUE)
cg_genelist=grep("^CAG", gene_list, value = TRUE)
ca_genelist=grep("^orf", gene_list, value = TRUE)
sc_genelist=grep("^Y", gene_list, value = TRUE)



## Sc GO annotate
##########################
go_annotations_sc <- gconvert(query = sc_genelist, 
                           organism = "scerevisiae", 
                           target = "GO", 
                           filter_na = TRUE)

# Create a named list where each element is a vector of GO terms for a gene
gene2GO_sc <- go_annotations_sc[,c("input","target")] %>%
  group_by(input) %>%
  summarize(target = list(target)) %>%
  deframe()

# Append any unsearched GO to the list
unsearched=sc_genelist[which(!(sc_genelist %in% go_annotations_sc$input))]
for(a in unsearched){
  path=paste0("D:/Grad study/lab/data/GO terms/",a,"_go_annotations.txt")
  x=read.table(path,sep="\t",skip=7,header=TRUE)[,c(1,3)]
  gene2GO_sc=append(gene2GO_sc,x %>%
           group_by(Gene.Format.Name) %>%
           summarize(target = list(Gene.Ontology.Term.ID)) %>%
           deframe())
}

# Save this as a data frame for future reference if needed
gene_GO_sc_df <- data.frame(
  Gene = rep(names(gene2GO_sc), sapply(gene2GO_sc, length)),
  GO = unlist(gene2GO_sc)
)
write.csv(gene_GO_sc_df,"D:/Grad study/lab/data/GO terms/annotated GO/gene_GO_sc.csv",
          row.names = FALSE)

##################################


## Cg GO annotate
#################################
go_annotations_cg <- gconvert(query = cg_genelist, 
                              organism = "cglabrata", 
                              target = "GO", 
                              filter_na = TRUE)

# Create a named list where each element is a vector of GO terms for a gene
gene2GO_cg <- go_annotations_cg[,c("input","target")] %>%
  group_by(input) %>%
  summarize(target = list(target)) %>%
  deframe()

# Append any unsearched GO 
unsearched=cg_genelist[which(!(cg_genelist %in% go_annotations_cg$input))]
cg_unsearched_go=read.csv("D:/Grad study/lab/data/GO terms/cg_rest.csv")
gene_GO_cg_unsearched <- data.frame(
  Gene = unlist(sapply(cg_unsearched_go$Gene.s.,function(x){strsplit(x," ")[[1]]})),
  GO = rep(cg_unsearched_go$GOID, 
           sapply(cg_unsearched_go$Gene.s.,function(x){
             length(strsplit(x," ")[[1]])
           }))
)
gene_GO_cg_unsearched$Gene[gene_GO_cg_unsearched$Gene=="DIG1"]="CAGL0L12782g"
gene_GO_cg_unsearched$Gene[gene_GO_cg_unsearched$Gene=="RME1"]="CAGL0K04257g"

gene_GO_cg_unsearched$GO=sapply(gene_GO_cg_unsearched$GO,function(x){
  paste0(c("GO:",rep(0,each=7-nchar(x)),x),collapse="")
})


# Save this as a data frame for future reference if needed
gene_GO_cg_df <- data.frame(
  Gene = rep(names(gene2GO_cg), sapply(gene2GO_cg, length)),
  GO = unlist(gene2GO_cg)
)
gene_GO_cg_df=rbind(gene_GO_cg_df,gene_GO_cg_unsearched)
write.csv(gene_GO_cg_df,"D:/Grad study/lab/data/GO terms/annotated GO/gene_GO_cg.csv",
          row.names = FALSE)
####################################

## Ca GO annotate
#################################
ca_go=read.csv("D:/Grad study/lab/data/GO terms/ca_rest.csv")

gene_GO_ca <- data.frame(
  Gene = unlist(sapply(ca_go$Gene.s.,function(x){strsplit(x," ")[[1]]})),
  GO = rep(ca_go$GOID, 
           sapply(ca_go$Gene.s.,function(x){
             length(strsplit(x," ")[[1]])
           }))
)

library(rvest)

get_orf_name <- function(gene) {
  tryCatch(
    {
      # get url
      url <- paste0("http://www.candidagenome.org/cgi-bin/locus.pl?locus=",gene,"&organism=C_albicans_SC5314")
    
    # read url
    webpage <- read_html(url)
      orf_name <- webpage %>%
      html_nodes(xpath = "//td[contains(., 'Assembly 19/21 Identifier')]/following-sibling::td") %>%
      html_text(trim = TRUE)
    return(orf_name[2])
    },
    error=function(cond) {
      message(gene)
      message("Here's the original error message:")
      message(conditionMessage(cond))
      # Choose a return value in case of error
      return(gene)
    }
  )
}
orf_results <- sapply(gene_GO_ca$Gene, get_orf_name)
# then search one by one the remaining genes.... and manually annotate
all_ca_correspond_orf=orf_results[which(grepl("orf",orf_results))]
ca_orf_id=all_ca_correspond_orf[!duplicated(all_ca_correspond_orf)]
gene_GO_ca$Gene=ca_orf_id[gene_GO_ca$Gene]


gene_GO_ca$GO=sapply(gene_GO_ca$GO,function(x){
  paste0(c("GO:",rep(0,each=7-nchar(x)),x),collapse="")
})

# Save this as a data frame for future reference if needed
write.csv(gene_GO_ca,"D:/Grad study/lab/data/GO terms/annotated GO/gene_GO_ca.csv",
          row.names = FALSE)

####################################

## Kl GO annotate: not done
####################################
kl_GO_raw=read.csv("D:/Grad study/lab/data/GO terms/annotate for Kl/E009_kl_GO_terms.csv")
kl_GO_raw=kl_GO_raw[grepl("KLL",kl_GO_raw$Gene.Names),c("Gene.Names","Gene.Ontology.IDs")]
kl_geneid=strsplit(kl_GO_raw$Gene.Names, split = "[ /;]")
kl_geneid=lapply(kl_geneid,function(x){x[grep("KLL",x)]})
gene=c()
go=c()
for(i in 1:length(kl_geneid)){
  
  go=append(go,rep(strsplit(kl_GO_raw$Gene.Ontology.IDs[i],"; ")[[1]],
            length(kl_geneid[[i]])))
  gene=append(gene,rep(kl_geneid[[i]],
                       each=length(strsplit(kl_GO_raw$Gene.Ontology.IDs[i],"; ")[[1]])))
}
cleaned_gene <- gsub("_", "", gene)
write.csv(data.frame(Gene=cleaned_gene,GO=go),"D:/Grad study/lab/data/GO terms/annotated GO/gene_GO_kl.csv",
          row.names = FALSE)


####################################

