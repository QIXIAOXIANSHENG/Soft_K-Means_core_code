gauss_bio=read.csv("D:/Grad study/lab/data/soft kmeans/gene selection/ygob_lbt001_cluster_result/topgo_result/without Kl/clustering scores/gaussian_bio_score.csv")
stat_score=read.csv("D:/Grad study/lab/data/soft kmeans/gene selection/ygob_lbt001_cluster_result/topgo_result/without Kl/clustering scores/stat_scores.csv")
db=stat_score$D.B

# db ratio
#######################################
db_r=db/max(db)
lambda=0.5
score_gauss_r_0.5=as.data.frame(as.matrix(gauss_bio[,-c(1,2)])*lambda+db_r)
score_gauss_r_0.5$k=gauss_bio$k
score_gauss_r_0.5$w=gauss_bio$w
score_gauss_r_0.5=score_gauss_r_0.5[,c(8,9,1:7)]

lambda=1
score_gauss_r_1=as.data.frame(as.matrix(gauss_bio[,-c(1,2)])*lambda+db_r)
score_gauss_r_1$k=gauss_bio$k
score_gauss_r_1$w=gauss_bio$w
score_gauss_r_1=score_gauss_r_1[,c(8,9,1:7)]


lambda=5
score_gauss_r=as.data.frame(as.matrix(gauss_bio[,-c(1,2)])*lambda+db_r)
score_gauss_r$k=gauss_bio$k
score_gauss_r$w=gauss_bio$w
score_gauss_r=score_gauss_r[,c(8,9,1:7)]


plot(x=1:310,y=score_gauss$score_self,type = "l",col="blue",ylim=c(0,35))
lines(x=1:310,y=score_gauss_r_1$score_self,col="red")
lines(x=1:310,y=score_gauss_r_0.5$score_self,col="black")

####################################

# db true value
####################################

lambda=0.5
score_gauss_v_0.5=as.data.frame(as.matrix(gauss_bio[,-c(1,2)])*lambda+db)
score_gauss_v_0.5$k=gauss_bio$k
score_gauss_v_0.5$w=gauss_bio$w
score_gauss_v_0.5=score_gauss_v_0.5[,c(8,9,1:7)]

lambda=1
score_gauss_v_1=as.data.frame(as.matrix(gauss_bio[,-c(1,2)])*lambda+db)
score_gauss_v_1$k=gauss_bio$k
score_gauss_v_1$w=gauss_bio$w
score_gauss_v_1=score_gauss_v_1[,c(8,9,1:7)]


lambda=5
score_gauss_v=as.data.frame(as.matrix(gauss_bio[,-c(1,2)])*lambda+db)
score_gauss_v$k=gauss_bio$k
score_gauss_v$w=gauss_bio$w
score_gauss_v=score_gauss_v[,c(8,9,1:7)]


plot(x=1:310,y=score_gauss_v_0.5$score_self,type = "l")
lines(x=1:310,y=score_gauss_v$score_self,col="blue")
lines(x=1:310,y=score_gauss_v_1$score_self,col="red")
