setwd("D:/source/RNASEQTest/RNASEQtest")

#load script for DE
source("./Rscript/DE.r")

#load script for normalization methods
source("./Rscript/NormMethods.r")


load("./dataset/HapMap.Rdata")

library(edgeR)

calCPM<-function(x,cond)
{
  nf <- calcNormFactors(x)
  lib.size<-colSums(x) * nf
  y <- t(log2(t(x + 0.5)/(lib.size + 1) * 1e+06))
  condA <- cond==cond[1]
  condB <- !condA
  yy<-pmax(rowMeans(y[,condA]),rowMeans(y[,condB]));
  names(yy)<- rownames(x)
  return(yy)
}
##calculate ordinary log2 foldchange
ordinaryFD <- function(x, condA, condB, log=TRUE)
{
  ma <- rowMeans(x[,condA])
  mb <- rowMeans(x[,condB])
  if(log) fc <- log2(mb+1)-log2(ma+1)
  else  fc <- mb-ma
  return (fc)
}

cond <- c(rep("condA",17),rep("condB",24))
##aFold
hapqt <- ResaFold(HapMap,cond)
ordqt <- ordinaryFD(qtotalnorm(HapMap),1:17,18:41)

##DESEQ2

hapds <- ResDESeq2(HapMap,cond)
ordds <- ordinaryFD(deseqnorm(HapMap),1:17,18:41)
##cpm
lcp <- calCPM(HapMap,cond)


cols <-rep("black",length(lcp))
cols[hapqt[,3]<.01] <- "red"
##plot
#rawFold DESEQ2

pdf(file="rc-hapmap-DESeq2.pdf",width=8,height=6)

plot(lcp,ordds,pch=16,ylim=c(-2.5,2.5),main="",ylab="",xlab="",axes=F,frame=T,col=cols)
axis(1,at=seq(-5,15,5),labels=rep("",5))
axis(2,at=seq(-3,3,1),labels=F)
abline(h=0.5,col="grey",lwd=2,lty=2)
abline(h=-0.5,col="grey",lwd=2,lty=2)
dev.off()

##shrinked lfc
pdf(file="sc-DESeq2.pdf",width=8,height=6)

plot(lcp,hapds[,1],pch=16,ylim=c(-2.5,2.5),main="",ylab="",xlab="",axes=F,frame=T,col=cols)

axis(1,at=seq(-5,15,5),labels=rep("",5))
axis(2,at=seq(-3,3,1),labels=F)
abline(h=0.5,col="grey",lwd=2,lty=2)
abline(h=-0.5,col="grey",lwd=2,lty=2)
dev.off()

#raw fc afold

pdf(file="rc-hapmap-afold.pdf",width=8,height=6)

plot(lcp,ordqt,pch=16,ylim=c(-2.5,2.5),main="",ylab="",xlab="",axes=F,frame=T,col=cols)

axis(1,at=seq(-5,15,5),labels=rep("",5))
axis(2,at=seq(-3,3,1),labels=F)
abline(h=0.5,col="grey",lwd=2,lty=2)
abline(h=-0.5,col="grey",lwd=2,lty=2)
dev.off()


pdf(file="sc-aFold.pdf",width=8,height=6)
#cols=rep("black",length(yy))
#cols[abs(res[["foldChange"]])>0.4]="red"
plot(lcp,hapqt[,1],pch=16,ylim=c(-2.5,2.5),main="",ylab="",xlab="",axes=F,frame=T,col=cols)

axis(1,at=seq(-5,15,5),labels=rep("",5))
axis(2,at=seq(-3,3,1),labels=F)
abline(h=0.5,col="grey",lwd=2,lty=2)
abline(h=-0.5,col="grey",lwd=2,lty=2)
dev.off()

###get DE for all methods on HapMap

cond <- c(rep("condA",17),rep("condB",24))
##aFold
hapqt <- ResaFold(HapMap,cond)

##DESEQ2
hapds <- ResDESeq2(HapMap,cond)

##voom
hapvm <- ResVoom(HapMap,cond)

##edgeR
haped <- ResedgeR(HapMap,cond)

##ABSSeq
hapab <- ResABSSeq(HapMap,cond)

##ROTS-total
x <- log2(totalnorm(HapMap)+1)
haprot <- ResROTS(x,cond)

##ROTS-qtotal
x <- log2(qtotalnorm(HapMap)+1)
haproq <- ResROTS(x,cond)


##baySeq
haprbs <- ResbaySeq(HapMap,cond)


