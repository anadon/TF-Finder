.libPaths("/var/www/html/cluster/R/x86_64-unknown-linux-gnu-library/3.0")
library(spls)

args<-commandArgs(TRUE)


TF<-read.csv(args[1],sep = "\t",header = F,stringsAsFactors = F)
# TF<-read.csv("/var/www/html/cluster/sample_data/TF-finder/allTF1640.txt",sep = "\t",header = F,stringsAsFactors = F)
oneGene <- TF[sample(1:dim(TF)[1], 1), ]
colnames(oneGene)<-NULL
write.table(oneGene,row.names = F,file = "predicted_pTFs.txt",sep = "\t",quote = FALSE)