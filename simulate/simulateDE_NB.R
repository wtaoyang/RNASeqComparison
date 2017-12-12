##moderate from deseq2
setwd("D:/source/RNASEQTest/RNASEQtest")
source("./Rscript/DE.r")

#load script for ROC 
source("./Rscript/ROC.r")
#load script for ROC 
source("./Rscript/NormMethods.r")

load("NB625res.Rdata")
load("NB2000res.Rdata")
load("NB1250res.Rdata")
load("NB4000res.Rdata")
load("NBR625res.Rdata")

###boxplot afold
pdf("NB_afold-2.pdf",width=12,height=8)

pd=cbind(NBres[["aFold"]][1:10,3],NBres[["aFoldt"]][1:10,3],NBres[["aFoldg"]][1:10,3],NBres[["aFoldq"]][1:10,3],
         NBRres[["aFold"]][1:10,3],NBRres[["aFoldt"]][1:10,3],NBRres[["aFoldg"]][1:10,3],NBRres[["aFoldq"]][1:10,3],
         NB2000res[["aFold"]][1:10,3],NB2000res[["aFoldt"]][1:10,3],NB2000res[["aFoldg"]][1:10,3],NB2000res[["aFoldq"]][1:10,3],
         NB1250res[["aFold"]][1:10,3],NB1250res[["aFoldt"]][1:10,3],NB1250res[["aFoldg"]][1:10,3],NB1250res[["aFoldq"]][1:10,3],
         NB4000res[["aFold"]][1:10,3],NB4000res[["aFoldt"]][1:10,3],NB4000res[["aFoldg"]][1:10,3],NB4000res[["aFoldq"]][1:10,3])

boxplot(pd,at=c(1:4,6:9,11:14,16:19,21:24),frame=T,xlab="",ylab="",main="",axes=F,xlim=c(0.5,24.5),cex.axis=0.7,ylim=c(0.6,1),col=2)
axis(2,at=seq(0.4,1.0,0.1),labels=rep("",length(seq(0.7,1.0,0.1))))
axis(1,at=1:4,labels=F)

abline(h=seq(0.4,1.0,.1),lty=2,col="gray")
abline(v=c(5,10,15,20),lty=2,col="gray")

dev.off()



fac=rep(c("aFold","DESeq2","edgeR","Voom","baySeq","ABSSeq","ROTS"),rep(10,7))
y=pd[,1:7]
#colnames(y)=c("ABSSeq","edgeR","DESeq","DESeq2","limmaVoom")
b<-data.frame(y=unlist(list(y)),factors=fac)
bartlett.test(y~factors,data=b)
m1<-aov(y~factors,data=b)
summary(m1)
TukeyHSD(m1)

y=pd[,8:14]
#colnames(y)=c("ABSSeq","edgeR","DESeq","baySeq","DESeq2","limmaVoom","limmaQN")
b<-data.frame(y=unlist(list(y)),factors=fac)
bartlett.test(y~factors,data=b)
m2<-aov(y~factors,data=b)
summary(m2)
TukeyHSD(m2)


y=pd[,15:21]
#colnames(y)=c("ABSSeq","edgeR","DESeq","baySeq","DESeq2","limmaVoom","limmaQN")
b<-data.frame(y=unlist(list(y)),factors=fac)
bartlett.test(y~factors,data=b)
m3<-aov(y~factors,data=b)
summary(m3)
TukeyHSD(m3)




