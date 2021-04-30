rare = read.table("E:\\学术研究\\西北大学\\microbiome\\GSMhindgutforgut\\data work\\N00306-ML1608041-西北大学-郭松涛-16S\\delivery\\Project_GST_ML1608041_16S_Result\\Result\\01_OTU\\otu_table_with_taxonomy.txt", header=T, sep="\t", row.names=1, comment.char="")#check.names=F,
#qiime1 的数据

# 提取genus，删除无用
genus=c()
delete=c()
count = 1
for(i in rare[, length(rare)])
{
  tmp = unlist(strsplit(as.character(i), split=";|__"))
  genus = c(genus, paste(tmp[12], sep="_"))#genus = c(genus, paste(tmp[10], tmp[12], sep="_"))
  if(is.na(tmp[10]) | tmp[10] == "" | is.na(tmp[12]) | tmp[12] == "")
  {
    delete = c(delete, count)
  }
  count = count + 1
}

rare$genus = genus
rare = rare[, -(ncol(rare)-1)]  # 删除注释
rare2 = rare[-delete, ]  # 删除NA和""
rare2 = rare2[order(rare2$genus, decreasing=F),]
plist = unique(rare2$genus)

# 合并行，相同门
rare3 = data.frame(apply(rare2[rare2$genus==plist[1], c(1:(ncol(rare2)-1))], 2, sum))
colnames(rare3)[1] = plist[1]
for(i in 2:length(plist))
{
  tmp = apply(rare2[rare2$genus==plist[i], c(1:(ncol(rare2)-1))], 2, sum)
  rare3 = cbind(rare3, tmp)
  colnames(rare3)[i] = plist[i]
}
rare3 = data.frame(t(rare3))

# 数据归一化
norm = rare3
sample_sum = apply(rare3, 2, sum)
for(i in 1:nrow(rare3))
{
  for(j in 1:ncol(rare3))
  {
    norm[i, j] = rare3[i, j]/sample_sum[j]
  }
}
apply(norm, 2, sum)  # 检测
write.csv(norm, file="genus_rare_norm.csv")

write.csv(rare3,file="genus_rare_orig.csv")

###############################################################
# 提取phylum，删除无用
phylum=c()
delete=c()
count = 1
for(i in rare[, length(rare)])
{
  tmp = unlist(strsplit(as.character(i), split=";|__"))
  phylum = c(phylum, tmp[4])
  if(is.na(tmp[4]) | tmp[4] == "")
  {
    delete = c(delete, count)
  }
  count = count + 1
}
rare$phylum = phylum
rare = rare[, -(ncol(rare)-1)]  # 删除注释
rare2 = rare[-delete, ]  # 删除NA和""
rare2 = rare2[order(rare2$phylum, decreasing=F),]
plist = unique(rare2$phylum)

# 合并行，相同门
rare3 = data.frame(apply(rare2[rare2$phylum==plist[1], c(1:(ncol(rare2)-1))], 2, sum))
colnames(rare3)[1] = plist[1]
for(i in 2:length(plist))
{
  tmp = apply(rare2[rare2$phylum==plist[i], c(1:(ncol(rare2)-1))], 2, sum)
  rare3 = cbind(rare3, tmp)
  colnames(rare3)[i] = plist[i]
}
rare3 = data.frame(t(rare3))
# 数据归一化
norm = rare3
sample_sum = apply(rare3, 2, sum)
for(i in 1:nrow(rare3))
{
  for(j in 1:ncol(rare3))
  {
    norm[i, j] = rare3[i, j]/sample_sum[j]
  }
}
apply(norm, 2, sum)  # 检测
write.csv(norm, file="phylum_rare_norm2.csv")
write.csv(rare3, file="phylum_rare_orig.csv")

#从otu中提取各个level的丰度数据。end
library(picante)
alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson') #   #Gini-Simpson 指数
  Pielou <- Shannon / log(Richness, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
  }
  result
}
IDcatal<-read.delim("E:\\学术研究\\西北大学\\microbiome\\GSMhindgutforgut\\清稿\\回复_ new version\\submit\\correct\\temp.txt",header=F)#map.txt
hIDcols<-IDcatal[grep("hindgut|Hindgut",IDcatal$V2),]$V1
fIDcols<-IDcatal[grep("foregut|Foregut",IDcatal$V2),]$V1
otug<-read.csv('E:\\学术研究\\西北大学\\microbiome\\GSMhindgutforgut\\data work\\郭松涛47sam-qiime2\\normalize\\tax_summary_a\\genus_headered.csv',row.names=1,stringsAsFactors=F,header=T,check.names=F)
otu<-read.csv('E:\\学术研究\\西北大学\\microbiome\\GSMhindgutforgut\\data work\\郭松涛47sam-qiime2\\normalize\\otu_taxa_table_headered.csv',row.names=1,stringsAsFactors=F,header=T,check.names=F)
otu<-read.delim('genus_rare_orig.csv',row.names=1,stringsAsFactors=F,header=T,check.names=F,sep=",")


#otuf<-otu[,which(names(otu)%in%paste("X",as.character(fIDcols),sep=""))]#otu[,paste("X",as.character(fIDcols),sep="")]
#otuh<-otu[,which(names(otu)%in%paste("X",as.character(hIDcols),sep=""))]
otuf<-otug[,which(names(otug)%in%fIDcols)]#otu[,paste("X",as.character(fIDcols),sep="")]
otuh<-otug[,which(names(otug)%in%hIDcols)]
otu<-t(otu)
otuf<-t(otuf)
otuh<-t(otuh)
alpha_all <- alpha(otu, base = 2)
alpha_f <- alpha(otuf, base = 2)
alpha_h <- alpha(otuh, base = 2)
write.csv(alpha_f, 'E:\\学术研究\\西北大学\\microbiome\\GSMhindgutforgut\\data work\\郭松涛47sam-qiime2\\normalize\\tax_summary_a\\asv_f_alpha.csv', quote = FALSE)
write.csv(alpha_h, 'E:\\学术研究\\西北大学\\microbiome\\GSMhindgutforgut\\data work\\郭松涛47sam-qiime2\\normalize\\tax_summary_a\\asv_h_alpha.csv', quote = FALSE)
write.csv(alpha_all,'E:\\学术研究\\西北大学\\microbiome\\GSMhindgutforgut\\data work\\郭松涛47sam-qiime2\\normalize\\tax_summary_a\\genus_all_alpha.csv',quote=FALSE)
alpha_f$group<-"foregut"
alpha_h$group<-"hindgut"
rb_alpha_f_h<-rbind(alpha_f,alpha_h)
par(mfrow = c(1, 2),bty="o" )
boxplot(Shannon~group,data=rb_alpha_f_h,col = c('dimgrey', 'red'),ylab="Shannon index(ASV)")
boxplot(alpha_h$Shannon)
#————————————————
#版权声明：本文为CSDN博主「马志远的生信笔记」的原创文章，遵循CC 4.0 BY-SA版权协议，转载请附上原文出处链接及本声明。
#原文链接：https://blog.csdn.net/weixin_42480153/article/details/108246368