.libPaths("/var/www/html/cluster/R/x86_64-unknown-linux-gnu-library/3.0")
library(spls)
cutoff = 0.05

# args<-commandArgs(TRUE)
args <- commandArgs(trailingOnly = TRUE)

TF<-read.csv(args[1],sep = "\t",header = F,stringsAsFactors = F)
# TF<-read.csv("/var/www/html/cluster/sample_data/TF-finder/allTF1640.txt",sep = "\t",header = F,stringsAsFactors = F)
tf_list<-as.character(TF$V1)
TF$V1<-NULL
TF<-as.data.frame(t(TF))
colnames(TF)<-tf_list
rownames(TF) <- NULL

# PW<-read.csv("/var/www/html/cluster/sample_data/TF-finder/target_lignin.txt",sep = "\t",header = F,stringsAsFactors = F)
PW<-read.csv(args[2],sep = "\t",header = F,stringsAsFactors = F)
pw_list<-as.character(PW$V1)
PW$V1<-NULL
PW<-as.data.frame(t(PW))
colnames(PW)<-pw_list
rownames(PW) <- NULL
cv <-cv.spls(TF,PW,eta = seq(0.6,0.9,0.1),K = c(5:10), 
             kappa=0.5, select="simpls", fit="kernelpls",scale.x=TRUE, scale.y=TRUE, plot.it=F)

eta = cv$eta.opt
K = cv$K.opt
f.ori <- spls(TF,PW,eta = eta,K = K)
ci.f <- ci.spls(f.ori,coverage = 0.95,plot.it = F,plot.var = F)
f.cor.mv <- correct.spls(ci.f,plot.it = F)
coef.f.corrected.mv <- abs(f.cor.mv)
sorted_score.cor.mv = data.frame(matrix(nrow = dim(TF)[2],ncol = 0))
pw_genes <-colnames(PW)

for (pw_i in 1:length(pw_genes)){
  temp<-data.frame(rownames(coef.f.corrected.mv),coef.f.corrected.mv[,pw_i])
  colnames(temp)<-c(paste(pw_genes[pw_i],"pwg",sep = "_"),paste(pw_genes[pw_i],"coeff",sep = "_"))
  temp<-temp[order(-temp[2]),]
  sorted_score.cor.mv <- cbind(sorted_score.cor.mv,temp)
  
}


Selected_TF<-data.frame(sort(table(unlist(sorted_score.cor.mv[1:10,seq(1,dim(sorted_score.cor.mv)[2],2)], use.names=FALSE)),decreasing = T))
#Selected_TF<-data.frame(Selected_TF[Selected_TF[,2]>2,])
Selected_TF<-data.frame(Selected_TF[Selected_TF[,1]>2,])
#Selected_TF <- Selected_TF$Var1
Selected_TF <- rownames(Selected_TF)
tf.exp<-TF[as.factor(Selected_TF)]
tf.exp<-t(tf.exp)
tf.exp<-as.data.frame(tf.exp)
tf.exp<-cbind(data.frame(Selected_TF),tf.exp)
colnames(tf.exp)<-NULL

write.table(tf.exp,row.names = F,file = "predicted_pTFs.txt",sep = "\t",quote = FALSE)
