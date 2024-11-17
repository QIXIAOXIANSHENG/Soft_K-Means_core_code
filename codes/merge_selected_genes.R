for(i in c("sc","ca","cg","kl")){
  assign(paste0(i,"_raw"),read.csv(paste0("D:\\Grad study\\lab\\data\\soft kmeans\\gene selection\\results_LBtest\\libnorm_",i,"_lbt_p_005.csv")))
  assign(i,
         as.matrix(get(paste0(i,"_raw"))[,seq(2,ncol(get(paste0(i,"_raw"))),2)])+as.matrix(get(paste0(i,"_raw"))[,seq(3,ncol(get(paste0(i,"_raw"))),2)]))
}

for(i in c("sc","ca","cg","kl")){
  assign(i,as.data.frame(get(i)))
}

## ca rowname: orf for ygob
ca_key=read.csv("D:/Grad study/lab/data/soft kmeans/orthogroup/Ca/Ca_gid_tid_orf.csv")
ca_key=ca_key[,-1]
ca$Geneid=ca_raw$X
names(ca_key)[2]="Geneid"
a=merge(ca_key,ca,by="Geneid")
ca=a[,c(3:16)]
names(ca)[1]="Geneid"
## kl rowname: KLLA0C12007g for ygob

new_id=c()
for(i in 1:nrow(kl)){
  new_id[i]=paste0(strsplit(kl_raw$X[i],"_")[[1]],collapse = "")
}
kl$Geneid=new_id

## Cg rownames: CAGL0A02904g for ygob
cg_key=read.csv("D:/Grad study/lab/data/soft kmeans/orthogroup/Cg/Cg_g_gid_p_n.csv")
cg_key=cg_key[,-1]
names(cg_key)[1]="Geneid"
cg$Geneid=cg_raw$X
a=merge(cg_key,cg,by="Geneid")
cg=a[,c(3:14)]
names(cg)[1]="Geneid"

## sc rownames: YDR362C for ygob
sc$Geneid=sc_raw$X


for(i in c("sc","ca","cg","kl")){
  print(colnames(get(i)))
}
# "ScB42"   "Sc30m2"  "Sc1h2"   "Sc1.5h2" "Sc2h2"   "Sc2.5h2" "Sc3h2"   "Sc3.5h2" "Sc4h2"   "Sc6h2"   "Sc8h2"  
# "CaB41"   "Ca15m1"  "Ca30m1"  "Ca45m1"  "Ca1h1"   "Ca1.5h1" "Ca2h1"   "Ca2.5h1" "Ca3h1"   "Ca3.5h1" "Ca4h1"   "Ca6h1"   "Ca8h1"  
# "CgB41"   "Cg15m1"  "Cg30m1"  "Cg45m1"  "Cg1h1"   "Cg1.5h1" "Cg2h1"   "Cg3h1"   "Cg4h1"   "Cg6h1"   "Cg8h1"  
# "KlB41"   "Kl15m3"  "Kl30m3"  "Kl45m3"  "Kl1h3"   "Kl1_5h3" "kl2h2"   "Kl2_5h1" "Kl3h3"   "Kl3_5h3" "Kl4h3"   "Kl6h3"   "Kl8h3"  

# B,30,1,1.5,2,3,4,6,8

sc=sc[,c(c(12,1,2,3,4,5,7,9,10,11))]
ca=ca[,c(c(0,1,3,5,6,7,9,11,12,13)+1)]
cg=cg[,c(c(0,1,3,5,6,7,8,9,10,11)+1)]
kl=kl[,c(c(14,1,3,5,6,7,9,11,12,13))]

colnames(sc)<-colnames(ca)<-colnames(cg)<-colnames(kl)<-c("Geneid","B","0.5","1","1.5","2","3","4","6","8")

all_4_species=rbind(sc, ca,cg,kl)
write.csv(all_4_species,"D:\\Grad study\\lab\\data\\soft kmeans\\gene selection\\ygob_libnorm_lbt005_all_4_species.csv",row.names = FALSE)
