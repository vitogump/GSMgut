#install.packages("fdrtool")
library("fdrtool")
#past wilcoxonpvalue col mannually
pathwaypvalue<-read.delim("E:\\学术研究\\西北大学\\microbiome\\GSMhindgutforgut\\清稿\\回复_ new version\\submit\\correct\\fdrtest.txt",header=T)
KOpvalue   <-  read.delim("E:\\学术研究\\西北大学\\microbiome\\GSMhindgutforgut\\清稿\\回复_ new version\\submit\\correct\\fdrtest.txt",header=T)
pwfdr = fdrtool(pathwaypvalue$wilcoxonpvalue, statistic="pvalue", plot =FALSE)
write.table(cbind(pwfdr$qval,pwfdr$lfdr),file="ssss",sep="\t",quote=F,row.names=F)

pwfdr = fdrtool(KOpvalue$wilcoxonpvalue, statistic="pvalue", plot =FALSE)
write.table(cbind(pwfdr$qval,pwfdr$lfdr),file="ssss",sep="\t",quote=F,row.names=F)
