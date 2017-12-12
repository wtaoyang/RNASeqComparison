setwd("D:/source/RNASEQTest/RNASEQtest")

#load script for DE
source("./Rscript/DE.r")

#load script for normalization methods
source("./Rscript/NormMethods.r")


load("./dataset/SEQC.Rdata")

##logCPM
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

##mix UHR and BHR
mixData=function(UHR,SEQC,normmethod)
{
  lfc <- abs(ordinaryFD(normmethod(SEQC),1:5,6:10)) ##get log FC
  slfc <- sample(c(lfc,-lfc),nrow(UHR)) #resample FC
  message(sum(slfc>0),"\t",sum(slfc<0))
  upl <- 2^pmax(0,slfc)
  dnl <- 2^(-pmin(0,slfc)) # ensure the total reads number are equally scaled
  dout <- cbind(UHR[,1:2]*upl,UHR[,3:4]*dnl)
  rout <- round(dout)
  dif <- dout -rout
  rout[dif>0.5] <- rout[dif>0.5]+1
  out <-rbind(as.matrix(UHR),as.matrix(rout))
  rownames(out) <- 1:nrow(out)
  return(out)
}

#resample
UHRi <- sample(1:5,4)
UHR <-SEQC[,UHRi]
ind <-rowSums(UHR)>0
UHR <-UHR[ind,]
BHR <-SEQC[rowSums(SEQC)>0,]
UHRqt <- mixData(UHR,BHR,qtotalnorm)
UHRds <- mixData(UHR,BHR,deseqnorm)

cond <- c("condA","condA","condB","condB")
##aFold
UHRaf <- ResaFold(UHR,cond)
UHRafm <- ResaFold(UHRqt,cond)
ordqt <- ordinaryFD(qtotalnorm(UHR),1:2,3:4)

##DESEQ2

UHRds <- ResDESeq2(UHR,cond)
UHRsdm <- ResDESeq2(UHRqt,cond)
ordds <- ordinaryFD(deseqnorm(UHR),1:2,3:4)
##cpm
lcp <- calCPM(UHR,cond)


##plot
#rawFold DESEQ2
pdf(file="rawfold-DESeq2.pdf",width=8,height=6)
plot(lcp,ordds,pch=16,ylim=c(-3,3),main="",ylab="",xlab="",axes=F,frame=T)
axis(1,at=seq(-5,15,5),labels=rep("",5))
axis(2,at=seq(-3,3,1),labels=F)
abline(h=0.5,col="grey",lwd=2,lty=2)
abline(h=-0.5,col="grey",lwd=2,lty=2)
dev.off()

##UHR
pdf(file="s1-DESeq2.pdf",width=8,height=6)

plot(lcp,UHRds[,1],pch=16,ylim=c(-1,1),main="",ylab="",xlab="",axes=F,frame=T)

axis(1,at=seq(-5,15,5),labels=rep("",5))
axis(2,at=seq(-1,1,0.5),labels=F)
abline(h=0.5,col="grey",lwd=2,lty=2)
abline(h=-0.5,col="grey",lwd=2,lty=2)
dev.off()

##UHR mixed with BHR
pdf(file="s2-DESeq2.pdf",width=8,height=6)

sx=lcp
sy=UHRsdm[1:nrow(UHR),1]
plot(sx,sy,pch=16,ylim=c(-3,3.2),main="",ylab="",xlab="",axes=F,frame=T)
axis(1,at=seq(-5,15,5),labels=rep("",5))
axis(2,at=seq(-3,3,1),labels=F)
abline(h=0.5,col="grey",lwd=2,lty=2)
abline(h=-0.5,col="grey",lwd=2,lty=2)
dev.off()

pdf(file="rawfold-afold.pdf",width=8,height=6)

plot(lcp,ordqt,pch=16,ylim=c(-3,3),main="",ylab="",xlab="",axes=F,frame=T)

axis(1,at=seq(-5,15,5),labels=rep("",5))
axis(2,at=seq(-3,3,1),labels=F)
abline(h=0.5,col="grey",lwd=2,lty=2)
abline(h=-0.5,col="grey",lwd=2,lty=2)
dev.off()


pdf(file="s1-aFold.pdf",width=8,height=6)

plot(lcp,UHRaf[,1],pch=16,ylim=c(-1,1),main="",ylab="",xlab="",axes=F,frame=T)

axis(1,at=seq(-5,15,5),labels=rep("",5))
axis(2,at=seq(-1,1,0.5),labels=F)
abline(h=0.5,col="grey",lwd=2,lty=2)
abline(h=-0.5,col="grey",lwd=2,lty=2)
dev.off()

pdf(file="s2-afold.pdf",width=8,height=6)
plot(lcp,UHRafm[1:nrow(UHR),1],pch=16,ylim=c(-1,1),main="",ylab="",xlab="",axes=F,frame=T)
axis(1,at=seq(-5,15,5),labels=rep("",5))
axis(2,at=seq(-1,1,0.5),labels=F)
abline(h=0.5,col="grey",lwd=2,lty=2)
abline(h=-0.5,col="grey",lwd=2,lty=2)
dev.off()




