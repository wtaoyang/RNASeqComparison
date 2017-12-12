setwd("D:/source/RNASEQTest/RNASEQtest")

#load script for normalization
source("./Rscript/NormMethods.r")

#load script for ROC 
source("./Rscript/ROC.r")


#load PrimePCR data to define true and false positives
load("./dataset/PrimePCR.Rdata")

#load ERCC data to calculate RMSD
load("./dataset/ERCCinfo.Rdata")
erccid <- rownames(ERCCinfo)

#gname in ensembl ID
qpcID <- rownames(PrimePCR)
#Tp
pset <- qpcID[abs(PrimePCR[,1])>0.5]
nset <- qpcID[abs(PrimePCR[,1])<0.2]

##calculate ordinary log2 foldchange
ordinaryFD <- function(x, condA, condB, log=TRUE)
{
  ma <- rowMeans(x[,condA])
  mb <- rowMeans(x[,condB])
  if(log) fc <- log2(mb+1)-log2(ma+1)
  else  fc <- mb-ma
  return (fc)
}
##calculate RMSD
RMSD <- function(x,y)
{
  return (sqrt(sum((x-y)^2)/length(x)))
}

#ABRF 
#load ABRF data set
load("./dataset/MAQC.Rdata")
xd <- MAQC
gname <- MAQCensid
gnamefilter <- gname[gname%in%qpcID]
condA <- 1:7
condB <- 8:14
xd <- xd[,c(8:14,1:7)]
##filter PrimePCR data

pset <- pset[pset %in% gname]
nset <- nset[nset %in% gname]
np <- length(pset)
nf <- length(nset)
statsCal <- function(tfc)
{
  ord <- order(abs(tfc),decreasing = T)
  roc<- roccurve(gname[ord],pset,nset)
  roc$erccRM <- RMSD(tfc[erccid],-ERCCinfo[erccid,5])
  roc$RM <- RMSD(tfc[gnamefilter],PrimePCR[gnamefilter,1])
  roc$se <- se(roc$auc,np,nf)
  
  return(roc)
}

##total
tmat <- totalnorm(xd)
tfc <- ordinaryFD(tmat,condA,condB)
names(tfc) <- gname
troc<- statsCal(tfc)

##uq
tmat <- upperqnorm(xd)
tfc <- ordinaryFD(tmat,condA,condB)
names(tfc) <- gname
uproc<- statsCal(tfc)

##qtotal
tmat <- qtotalnorm(xd)
tfc <- ordinaryFD(tmat,condA,condB)
names(tfc) <- gname
qtroc<- statsCal(tfc)
##tmm

tmat <- tmmnorm(xd)
tfc <- ordinaryFD(tmat,condA,condB)
names(tfc) <- gname
tmroc<- statsCal(tfc)

##DESEQ
tmat <- deseqnorm(xd)
tfc <- ordinaryFD(tmat,condA,condB)
names(tfc) <- gname
dsroc<- statsCal(tfc)
##cqn

tmat <- cqnnorm(xd,MAQCgcl[,1],MAQCgcl[,2])
tfc <- ordinaryFD(tmat,condA,condB,F)
names(tfc) <- gname
cqroc<- statsCal(tfc)

##MedpgQ2
tmat <- MedpgQ2(xd)
tfc <- ordinaryFD(tmat,condA,condB)
names(tfc) <- gname
medroc<- statsCal(tfc)

##uqpgQ2
tmat <- UQpgQ2(xd)
tfc <- ordinaryFD(tmat,condA,condB)
names(tfc) <- gname
uqroc<- statsCal(tfc)

###z-score of AUC compared to qtotal
Ztot <- (qtroc$auc-troc$auc)/sqrt(qtroc$se^2+troc$se^2)
Ztotmm <- (qtroc$auc-tmroc$auc)/sqrt(qtroc$se^2+tmroc$se^2)
Ztods <- (qtroc$auc-dsroc$auc)/sqrt(qtroc$se^2+dsroc$se^2)
Ztoup <- (qtroc$auc-uproc$auc)/sqrt(qtroc$se^2+uproc$se^2)
Ztomed <- (qtroc$auc-medroc$auc)/sqrt(qtroc$se^2+medroc$se^2)
Ztouq <- (qtroc$auc-uqroc$auc)/sqrt(qtroc$se^2+uqroc$se^2)
Ztocq <- (qtroc$auc-cqroc$auc)/sqrt(qtroc$se^2+cqroc$se^2)

zscore <- c(Ztotmm,Ztot,Ztoup,Ztods,Ztocq,Ztomed,Ztouq)
zp <- pnorm(zscore,lower.tail = F)*length(zscore)


aucs <- c(qtroc$auc,tmroc$auc,troc$auc,uproc$auc,dsroc$auc,cqroc$auc,medroc$auc,uqroc$auc)



##plot
pdf("MAQCAUC_norm.pdf",width=6,height=6)

par(mar=c(0,0,0,0)+5)
#
plot(qtroc,type="l",lwd=3,main="",xlab="",ylab="",col=2,axes=F,frame=T,xlim=c(0,1),ylim=c(0,1))
axis(1,at=seq(0,1,length=5),labels=FALSE)
axis(2,at=seq(0,1,length=5),labels=FALSE)
lines(tmroc,lwd=3,col=3)
lines(troc,lwd=3,col=4)
lines(uproc,lwd=3,col=5) 
lines(dsroc,lwd=3,col=6)
lines(cqroc,lwd=3,col=7)
lines(medroc,lwd=3,col=8)
lines(uqroc,lwd=3,col=1)
lines(qtroc,lwd=3,col=2)
lines(seq(0,1,0.05),seq(0,1,0.05),lwd=3,col="gray",lty=2)

dev.off()



pdf("MAQCAUC_RM.pdf",width=6,height=6)

par(mar=c(0,0,0,0)+5)
#
dat <- c(qtroc$RM,tmroc$RM,troc$RM,uproc$RM,
         dsroc$RM,cqroc$RM,medroc$RM,uqroc$RM)
barplot(dat,main="",xlab="",ylab="",col=c(2:8,1),axes=F,xlim=c(0,9.5),ylim=c(0,6))
axis(1,at=seq(0.7,9.3,1.2),labels=FALSE)
axis(2,at=seq(0,6,length=4),labels=FALSE)

dev.off()

