setwd("D:/source/RNASEQTest/RNASEQtest")

#load script for DE
source("./Rscript/DE.r")

#load script for ROC 
source("./Rscript/ROC.r")


#load PrimePCR data to define true and false positives
load("./dataset/PrimePCR.Rdata")

#load ERCC data to calculate RMSD
load("./dataset/ERCCinfo.Rdata")
erccid <- rownames(ERCCinfo)

#gname gene name
qpcID <- rownames(PrimePCR)
#Tp
pset <- qpcID[abs(PrimePCR[,1])>0.5]
nset <- qpcID[abs(PrimePCR[,1])<0.2]
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
##calculate RMSD
RMSD <- function(x,y)
{
  return (sqrt(sum((x-y)^2)/length(x)))
}


#ABRFAB 
#load ABRF data set
load("./dataset/ABRFAB.Rdata")
gname <- ABRFABensid
gnamefilter <- gname[gname%in%qpcID]
xd <- ABRFAB
condA <- 1:9
condB <- 10:18
cond <- factor(c(rep("condA",9),rep("condB",9)))

sfq <- qtotalNormalized(ABRFAB)
##filter PrimePCR data

pset <- pset[pset %in% gname]
nset <- nset[nset %in% gname]
np <- length(pset)
nf <- length(nset)
##calculate statistics
# fold, pvalue, adjusted pvalue
statsCal <- function(tfc,pv,adpv)
{
  ord <-c()
  if(missing(pv))
  {
    ord <- order(abs(tfc),decreasing = T)
  }else
  {
    ord <- order(pv,-abs(tfc),decreasing = F)
  }
  roc<- roccurve(gname[ord],pset,nset)
  roc$erccRM <- RMSD(tfc[erccid],-ERCCinfo[erccid,5])
  roc$RM <- RMSD(tfc[gnamefilter],PrimePCR[gnamefilter,1])
  roc$se <- se(roc$auc,np,nf)
  
  aset<- gname[adpv<0.05]
  x <- sum(aset%in%pset)#/np
  y <- sum(aset%in%nset)#/(#fncurve(gname[ord],pset,nset)
  roc$efdr <- y/(x+y)
  roc$sens <- x/np
  return(roc)
}

##aFold
tafold <- read.table("tafold-ABRF.txt",header=T,row.names=1,sep="\t")#ResaFold(xd,cond)
tfc <- tafold[,1]
names(tfc) <- gname
aroc <- statsCal(tfc,adpv=tafold[,3])
write.table(tafold,"tafold-ABRF.txt",quote=F,sep="\t")


##aFold TMM
tafoldtmm <- ResaFold(xd,cond,normM="TMM")
tfc <- tafoldtmm[,1]
names(tfc) <- gname
aroct <- statsCal(tfc,adpv=tafold[,3])

write.table(tafoldtmm,"tafold-tmm-ABRF.txt",quote=F,sep="\t")

##aFold geometric
tafoldgeo <- ResaFold(xd,cond,normM="geometric")
tfc <- tafoldgeo[,1]
names(tfc) <- gname
arocg <- statsCal(tfc,adpv=tafoldgeo[,3])

write.table(tafoldgeo,"tafold-geo-ABRF.txt",quote=F,sep="\t")

##aFold quartile
tafoldqua <- ResaFold(xd,cond,normM="quartile")
tfc <- tafoldqua[,1]
names(tfc) <- gname
arocq <- statsCal(tfc,adpv=tafoldqua[,3])

write.table(tafoldqua,"tafold-qua-ABRF.txt",quote=F,sep="\t")

##DESEQ2
tdeseq2 <- read.table("tdeseq2-ABRF.txt",header=T,row.names=1,sep="\t")#ResDESeq2(xd,cond)
tfc <- tdeseq2[,1]
pv <- tdeseq2[,2]
names(tfc) <- gname
dsroc <- statsCal(tfc,pv,adpv=tdeseq2[,3])
write.table(tdeseq2,"tdeseq2-ABRF.txt",quote=F,sep="\t")


##DESEQ2 qtotal
tdeseq2q <- ResDESeq2(xd,cond,sf=sfq)
tfc <- tdeseq2q[,1]
pv <- tdeseq2q[,2]
names(tfc) <- gname
dsrocq <- statsCal(tfc,pv,adpv=tdeseq2q[,3])
write.table(tdeseq2q,"tdeseq2-q-ABRF.txt",quote=F,sep="\t")
##edgeR
tedgeR <- read.table("tedgeR-ABRF.txt",header=T,row.names=1,sep="\t")#ResedgeR(xd,cond)
tedgeR <- tedgeR[rownames(xd),]
tfc <- tedgeR[,1]
pv <- tedgeR[,2]
names(tfc) <- gname
edroc <- statsCal(tfc,pv,adpv=tedgeR[,3])
write.table(tedgeR,"tedgeR-ABRF.txt",quote=F,sep="\t")
##edgeR qtotal
tedgeRq <- ResedgeR(xd,cond,sf=sfq)
tedgeRq <- tedgeRq[rownames(xd),]
tfc <- tedgeRq[,1]
pv <- tedgeRq[,2]
names(tfc) <- gname
edrocq <- statsCal(tfc,pv,adpv=tedgeRq[,3])
write.table(tedgeRq,"tedgeR-q-ABRF.txt",quote=F,sep="\t")
##voom

tVoom <- read.table("tVoom-ABRF.txt",header=T,row.names=1,sep="\t")#ResVoom(xd,cond)
tVoom <- tVoom[rownames(xd),]
tfc <- tVoom[,1]
names(tfc) <- gname
pv <- tVoom[,2]
vmroc <- statsCal(tfc,pv,adpv=tVoom[,3])
write.table(tVoom,"tVoom-ABRF.txt",quote=F,sep="\t")

##voom qtotal

tVoomq <- ResVoom(xd,cond,sf=sfq)
tVoomq <- tVoomq[rownames(xd),]
tfc <- tVoomq[,1]
names(tfc) <- gname
pv <- tVoomq[,2]
vmrocq <- statsCal(tfc,pv,adpv=tVoomq[,3])
write.table(tVoomq,"tVoom-q-ABRF.txt",quote=F,sep="\t")


##baySeq
#tbaySeq <- ResbaySeq(xd,cond,cl=8)
tbaySeq <- read.table("ABRFAB_res.txt",header=T,row.names=1,sep="\t")
tbaySeq <- tbaySeq[rownames(xd),]
tfc <- ordinaryFD(upperqnorm(xd),condA,condB)
names(tfc) <- gname
pv <- tbaySeq[,2]
bsroc <- statsCal(tfc,pv,adpv=pv)

##baySeq qtotal
#tbaySeqq <- ResbaySeq(xd,cond,cl=8)
tbaySeqq <- read.table("ABRFAB_res_qt.txt",header=T,row.names=1,sep="\t")
tbaySeqq <- tbaySeqq[rownames(xd),]
tfc <- ordinaryFD(qtotalnorm(xd),condA,condB)
names(tfc) <- gname
pv <- tbaySeqq[,2]
bsrocq <- statsCal(tfc,pv,adpv=pv)
##ABSSeq

tABSSeq <- read.table("tABSSeq-ABRF.txt",header=T,row.names=1,sep="\t")#ResABSSeq(xd,cond)
tfc <- tABSSeq[,1]
names(tfc) <- gname
pv <- tABSSeq[,2]
abroc <- statsCal(tfc,pv,adpv=tABSSeq[,3])
write.table(tABSSeq,"tABSSeq-ABRF.txt",quote=F,sep="\t")

##ROTS-total
source("./Rscript/NormMethods.r")
x <- log2(totalnorm(xd)+1)
tROTSt <- read.table("tROTSt-ABRF.txt",header=T,row.names=1,sep="\t")#ResROTS(x,cond)
tROTSt <- tROTSt[rownames(xd),]
tfc <- -tROTSt[,1]
names(tfc) <- gname
pv <- tROTSt[,2]
rotroc <- statsCal(tfc,pv,adpv=tROTSt[,3])
write.table(tROTSt,"tROTSt-ABRF.txt",quote=F,sep="\t")
##ROTS-qtotal
x <- log2(qtotalnorm(xd)+1)
qROTSt <- read.table("qROTSt-ABRF.txt",header=T,row.names=1,sep="\t")#ResROTS(x,cond)
qROTSt <- qROTSt[rownames(xd),]
tfc <- -qROTSt[,1]
names(tfc) <- gname
pv <- qROTSt[,2]
roqroc <- statsCal(tfc,pv,adpv=qROTSt[,3])

write.table(qROTSt,"qROTSt-ABRF.txt",quote=F,sep="\t")
###z-score of AUC compared to afold

Ztods <- abs(aroc$auc-dsroc$auc)/sqrt(aroc$se^2+dsroc$se^2)
Ztodsq <- abs(aroc$auc-dsrocq$auc)/sqrt(aroc$se^2+dsrocq$se^2)
Ztodsi <- abs(dsrocq$auc-dsroc$auc)/sqrt(dsroc$se^2+dsrocq$se^2)
Ztoed <- abs(aroc$auc-edroc$auc)/sqrt(aroc$se^2+edroc$se^2)
Ztoedq <- abs(aroc$auc-edrocq$auc)/sqrt(aroc$se^2+edrocq$se^2)
Ztoedi <- abs(edrocq$auc-edroc$auc)/sqrt(edrocq$se^2+edroc$se^2)
Ztovm <- abs(aroc$auc-vmroc$auc)/sqrt(aroc$se^2+vmroc$se^2)
Ztovmq <- abs(aroc$auc-vmrocq$auc)/sqrt(aroc$se^2+vmrocq$se^2)
Ztovmi <- abs(vmrocq$auc-vmroc$auc)/sqrt(vmrocq$se^2+vmroc$se^2)
Ztobs <- abs(aroc$auc-bsroc$auc)/sqrt(aroc$se^2+bsroc$se^2)
Ztobsq <- abs(aroc$auc-bsrocq$auc)/sqrt(aroc$se^2+bsrocq$se^2)
Ztobsi <- abs(bsrocq$auc-bsroc$auc)/sqrt(bsrocq$se^2+bsroc$se^2)
Ztoab <- abs(aroc$auc-abroc$auc)/sqrt(aroc$se^2+abroc$se^2)
Ztoro <- abs(aroc$auc-rotroc$auc)/sqrt(aroc$se^2+rotroc$se^2)
Ztoroq <- abs(aroc$auc-roqroc$auc)/sqrt(aroc$se^2+roqroc$se^2)
Ztoroi <- abs(roqroc$auc-rotroc$auc)/sqrt(roqroc$se^2+rotroc$se^2)

zscore <- c(Ztods,Ztodsq,Ztoed,Ztoedq,Ztovm,Ztovmq,Ztobs,Ztobsq,Ztoab,Ztoro
            ,Ztoroq,Ztodsi,Ztoedi,Ztovmi,Ztobsi,Ztoroi)
zp <- pnorm(zscore,lower.tail = F)*length(zscore)
write.table(cbind(zscore,zp),"ABRF_DE_all.txt",quote=F,sep="\t",row.names=F,col.names = F)


aucs <- c(aroc$auc,dsroc$auc,dsrocq$auc,edroc$auc,edrocq$auc,vmroc$auc,vmrocq$auc,
          bsroc$auc,bsrocq$auc,abroc$auc,rotroc$auc,roqroc$auc)

write.table(aucs,"ABRF_DE_all_AUC.txt",quote=F,sep="\t",row.names=F,col.names = F)


##afold switch

Ztodsg <- abs(arocg$auc-dsroc$auc)/sqrt(arocg$se^2+dsroc$se^2)

Ztoedt <- abs(aroct$auc-edroc$auc)/sqrt(aroct$se^2+edroc$se^2)

Ztovmt <- abs(aroct$auc-vmroc$auc)/sqrt(aroct$se^2+vmroc$se^2)

Ztobsqua <- abs(arocq$auc-bsroc$auc)/sqrt(arocq$se^2+bsroc$se^2)


zscore <- c(Ztodsg,Ztoedt,Ztovmt,Ztobsqua)
zp <- pnorm(zscore,lower.tail = F)*length(zscore)
write.table(cbind(zscore,zp),"ABRF_DE_afold_switch.txt",quote=F,sep="\t",row.names=F,col.names = F)


aucs <- c(aroct$auc,arocq$auc,arocg$auc)

write.table(aucs,"ABRF_DE_all_AUC_afold_swtich.txt",quote=F,sep="\t",row.names=F,col.names = F)

pdf("ABRFAUC_de.pdf",width=6,height=6)

par(mar=c(0,0,0,0)+5)
#
plot(aroc,type="l",lwd=3,main="",xlab="",ylab="",col=2,axes=F,frame=T,xlim=c(0,1),ylim=c(0,1))
axis(1,at=seq(0,1,length=5),labels=FALSE)
axis(2,at=seq(0,1,length=5),labels=FALSE)
lines(dsroc,lwd=3,col=3)
lines(dsrocq,lwd=3,col=3,lty=2)
lines(edroc,lwd=3,col=4) 
lines(edrocq,lwd=3,col=4,lty=2)
lines(vmroc,lwd=3,col=5)
lines(vmrocq,lwd=3,col=5,lty=2)
lines(bsroc,lwd=3,col=6)
lines(bsrocq,lwd=3,col=6,lty=2)
lines(abroc,lwd=3,col=7)
lines(rotroc,lwd=3,col=8)
lines(roqroc,lwd=3,col=8,lty=2)
lines(seq(0,1,0.05),seq(0,1,0.05),lwd=3,col="gray",lty=2)

lines(c(0.65,0.7),rep(0.52,2),lwd=3,col=2);#text(0.80,0.3,"aFold");
lines(c(0.65,0.7),rep(0.45,2),lwd=3,col=3);#text(0.79,0.25,"DESeq2");
lines(c(0.65,0.7),rep(0.38,2),lwd=3,col=4);#text(0.79,0.2,"edgeR");
lines(c(0.65,0.7),rep(0.31,2),lwd=3,col=5);#text(0.80,0.15,"Voom");
lines(c(0.65,0.7),rep(0.24,2),lwd=3,col=6);#text(0.79,0.10,"baySeq");
lines(c(0.65,0.7),rep(0.17,2),lwd=3,col=7);#text(0.79,0.05,"ABSSeq");
lines(c(0.65,0.7),rep(0.10,2),lwd=3,col=8);#text(0.79,0.35,"ROTS");

dev.off()

pdf("ABRFfn_de.pdf",width=6,height=6)

par(mar=c(0,0,0,0)+5)
#
plot(aroc$efdr,aroc$sens,lwd=3,main="",xlab="",ylab="",col=2,axes=F,frame=T,xlim=c(0,0.1),ylim=c(0.7,1),pch=16,cex=1.5)
axis(1,at=seq(0,.1,length=5),labels=FALSE)
axis(2,at=seq(0.7,1,length=4),labels=FALSE)
points(dsroc$efdr,dsroc$sens,col=3,pch=16,cex=1.5)
points(dsrocq$efdr,dsrocq$sens,col=3,pch=1,cex=1.5)
points(edroc$efdr,edroc$sens,col=4,pch=16,cex=1.5) 
points(edrocq$efdr,edrocq$sens,col=4,pch=1,cex=1.5)
points(vmroc$efdr,vmroc$sens,col=5,pch=16,cex=1.5) 
points(vmrocq$efdr,vmrocq$sens,col=5,pch=1,cex=1.5) 
points(bsroc$efdr,bsroc$sens,col=6,pch=16,cex=1.5) 
points(bsrocq$efdr,bsrocq$sens,col=6,pch=1,cex=1.5) 
points(abroc$efdr,abroc$sens,col=7,pch=16,cex=1.5) 
points(rotroc$efdr,rotroc$sens,col=8,pch=16,cex=1.5) 
points(roqroc$efdr,roqroc$sens,col=8,pch=1,cex=1.5)
abline(v=0.05,col="grey",lwd=2,lty=2)
abline(h=0.9,col="grey",lwd=2,lty=2)
#lines(c(0.65,0.7),rep(0.52,2),lwd=3,col=2);#text(0.80,0.3,"aFold");
#lines(c(0.65,0.7),rep(0.45,2),lwd=3,col=3);#text(0.79,0.25,"DESeq2");
#lines(c(0.65,0.7),rep(0.38,2),lwd=3,col=4);#text(0.79,0.2,"edgeR");
#lines(c(0.65,0.7),rep(0.31,2),lwd=3,col=5);#text(0.80,0.15,"Voom");
#lines(c(0.65,0.7),rep(0.24,2),lwd=3,col=6);#text(0.79,0.10,"baySeq");
#lines(c(0.65,0.7),rep(0.17,2),lwd=3,col=7);#text(0.79,0.05,"ABSSeq");
#lines(c(0.65,0.7),rep(0.10,2),lwd=3,col=8);#text(0.79,0.35,"ROTS");
#lines(c(0.65,0.7),rep(0.03,2),lwd=3,col=1);#text(0.79,0.35,"UppgQ2");
dev.off()


pdf("ABRFAUC_de_afold.pdf",width=6,height=6)

par(mar=c(0,0,0,0)+5)
#
plot(aroc,type="l",lwd=3,main="",xlab="",ylab="",col=2,axes=F,frame=T,xlim=c(0,1),ylim=c(0,1))
axis(1,at=seq(0,1,length=5),labels=FALSE)
axis(2,at=seq(0,1,length=5),labels=FALSE)
lines(aroct,lwd=3,col=2,lty=2)
lines(arocg,lwd=3,col=2,lty=3)
lines(arocq,lwd=3,col=2,lty=4)
lines(dsroc,lwd=3,col=3)
lines(edroc,lwd=3,col=4) 
lines(vmroc,lwd=3,col=5)
lines(bsroc,lwd=3,col=6)
#lines(abroc,lwd=3,col=7)

lines(seq(0,1,0.05),seq(0,1,0.05),lwd=3,col="gray",lty=2)

lines(c(0.65,0.7),rep(0.52,2),lwd=3,col=2);#text(0.80,0.3,"aFold");
lines(c(0.65,0.7),rep(0.45,2),lwd=3,col=2,lty=2);#text(0.79,0.25,"aFold-tmm");
lines(c(0.65,0.7),rep(0.38,2),lwd=3,col=2,lty=3);#text(0.79,0.2,"aFold-geometric");
lines(c(0.65,0.7),rep(0.31,2),lwd=3,col=2,lty=4);#text(0.80,0.15,"aFold-quartile");
lines(c(0.65,0.7),rep(0.24,2),lwd=3,col=3);#text(0.79,0.10,"DESEQ2");
lines(c(0.65,0.7),rep(0.17,2),lwd=3,col=4);#text(0.79,0.05,"edgeR");
lines(c(0.65,0.7),rep(0.10,2),lwd=3,col=5);#text(0.79,0.35,"Voom");
lines(c(0.65,0.7),rep(0.03,2),lwd=3,col=6);#text(0.79,0.35,"BaySeq");
dev.off()

pdf("ABRFAUC_ERCC_de.pdf",width=6,height=6)

par(mar=c(0,0,0,0)+5)
#
dat <- c(bsrocq$erccRM,aroc$erccRM,dsrocq$erccRM)
cols <- c("grey","black","black")
#dens <- c(1000,20,1000,20,50)
barplot(dat,main="",xlab="",ylab="",col=cols,axes=F,xlim=c(0,4),ylim=c(0,1))
axis(1,at=seq(0.7,3.3,1.2),labels=FALSE)
axis(2,at=seq(0,1,length=5),labels=FALSE)

dev.off()
pdf("ABRFAUC_RM_de.pdf",width=6,height=6)

par(mar=c(0,0,0,0)+5)
#
dat <- c(bsrocq$RM,aroc$RM,dsrocq$RM)
cols <- c("grey","black","black")

barplot(dat,main="",xlab="",ylab="",col=cols,axes=F,xlim=c(0,4),ylim=c(0,8))
axis(1,at=seq(0.7,3.3,1.2),labels=FALSE)
axis(2,at=seq(0,8,length=5),labels=FALSE)

dev.off()

pdf("ABRFAUC_RM_de_all.pdf",width=8,height=6)

par(mar=c(0,0,0,0)+5)
#
dat <- c(aroc$RM,dsroc$RM,edroc$RM,vmroc$RM,bsroc$RM,abroc$RM,rotroc$RM,roqroc$RM)
cols <- c(2:8,1)

barplot(dat,main="",xlab="",ylab="",col=cols,axes=F,xlim=c(0,9),ylim=c(0,8.5))
axis(1,at=seq(0.7,9.3,1.2),labels=FALSE)
axis(2,at=seq(0,8.5,length=5),labels=FALSE)

dev.off()

pdf("ABRFAUC_ercc_de_all.pdf",width=8,height=6)

par(mar=c(0,0,0,0)+5)
#
dat <- c(aroc$erccRM,dsroc$erccRM,edroc$erccRM,vmroc$erccRM,bsroc$erccRM,abroc$erccRM,rotroc$erccRM,roqroc$erccRM)
cols <- c(2:8,1)

barplot(dat,main="",xlab="",ylab="",col=cols,axes=F,xlim=c(0,9),ylim=c(0,1))
axis(1,at=seq(0.7,9.3,1.2),labels=FALSE)
axis(2,at=seq(0,1,length=4),labels=FALSE)

dev.off()
##get plot for FC from PrimePCR and ABRF
lcp <- calCPM(xd,cond)
names(lcp) <- gname
lcp <- lcp[gnamefilter]
cols<- rep("black",length(lcp))
cols[lcp<1] <- "red"
##fold change correlation
afd <- tafold[,1]
names(afd) <- gname
dsd <- tdeseq2q[,1]
names(dsd) <- gname
rwd <- ordinaryFD(qtotalnorm(xd),condA,condB)
names(rwd) <- gname


#cor(x[lcp>1],y[lcp>1])
pdf(file="ABRF-af-Prime.pdf",width=8,height=6)
x <-afd[gnamefilter]
y<-PrimePCR[gnamefilter,1]

plot(x,y,pch=16,ylim=c(-45,45),main="",ylab="",xlab="",axes=F,frame=T,col=cols)

axis(1,at=seq(-6,6,2),labels=F)
axis(2,at=seq(-40,40,20),labels=F)
dev.off()
pdf(file="ABRF-raw-Prime.pdf",width=8,height=6)
x <-rwd[gnamefilter]
y<-PrimePCR[gnamefilter,1]

plot(x,y,pch=16,ylim=c(-45,45),main="",ylab="",xlab="",axes=F,frame=T,col=cols)

axis(1,at=seq(-10,15,5),labels=F)
axis(2,at=seq(-40,40,20),labels=F)
dev.off()

pdf(file="ABRF-ds-Prime.pdf",width=8,height=6)
x <-dsd[gnamefilter]
y<-PrimePCR[gnamefilter,1]

plot(x,y,pch=16,ylim=c(-45,45),main="",ylab="",xlab="",axes=F,frame=T,col=cols)

axis(1,at=seq(-10,15,5),labels=F)
axis(2,at=seq(-40,40,20),labels=F)
dev.off()


