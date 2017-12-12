##moderate from deseq2
setwd("D:/source/RNASEQTest/RNASEQtest")
source("./Rscript/DE.r")

#load script for ROC 
source("./Rscript/ROC.r")
#load script for ROC 
source("./Rscript/NormMethods.r")
statsCal <- function(tfc,pv,adpv,beta,gname)
{
  idx<-adpv<0.1
  ord <-c()
  if(missing(pv))
  {
    ord <- order(abs(tfc),decreasing = T)
  }else
  {
    ord <- order(pv,-abs(tfc),decreasing = F)
  }
  roc <-roccurve(gname[ord],pset,nset)
  sens <- mean(idx[abs(beta)>0])
  idx <- which(adpv<0.1)
  efdr <- ifelse(sum(idx) == 0, 0, mean((beta == 0)[idx]))
  return(c(sens,efdr,roc$auc))
}


afres <- matrix(c(0),nrow=30,ncol=3)
afrest <- matrix(c(0),nrow=30,ncol=3)
afresg <- matrix(c(0),nrow=30,ncol=3)
afresq <- matrix(c(0),nrow=30,ncol=3)
dsres <- matrix(c(0),nrow=30,ncol=3)
dsresq <- matrix(c(0),nrow=30,ncol=3)
edres <- matrix(c(0),nrow=30,ncol=3)
edresq <- matrix(c(0),nrow=30,ncol=3)
vmres <- matrix(c(0),nrow=30,ncol=3)
vmresq <- matrix(c(0),nrow=30,ncol=3)
bsres <- matrix(c(0),nrow=30,ncol=3)
bsresq <- matrix(c(0),nrow=30,ncol=3)
abres <- matrix(c(0),nrow=30,ncol=3)
rtres <- matrix(c(0),nrow=30,ncol=3)
rtresq <- matrix(c(0),nrow=30,ncol=3)
basn="CT_2_"
path="./dataset/SimStudyCT_B_625_625/2samplespergroup"
cond<-c(rep("condA",2),rep("condB",2))
for(i in 1:10)
{
  filen<-paste(path,"/",basn,i,".rds",sep="")
  dat <-readRDS(filen)
  mat <-dat[[1]]
  gname <- rownames(mat)
  pset <- as.character(dat[[3]])
  indx <- gname%in%pset
  nset <- gname[!indx]
  beta <- rep(0,nrow(mat))
  beta[indx] <- log(1.5)
  sfs <- qtotalNormalized(mat)
  tres <- ResaFold(mat,cond)
  afres[i,] <-statsCal(tres[,1],adpv=tres[,3],beta=beta,gname=gname)
  #tres <- ResaFold(mat,cond,normM = "TMM")
  #afrest[i,] <-statsCal(tres[,1],adpv=tres[,3],beta=beta,gname=gname)
  #tres <- ResaFold(mat,cond,normM = "geometric")
  #afresg[i,] <-statsCal(tres[,1],adpv=tres[,3],beta=beta,gname=gname)
  #tres <- ResaFold(mat,cond,normM = "quartile")
  #afresq[i,] <-statsCal(tres[,1],adpv=tres[,3],beta=beta,gname=gname)
  #tres <- ResDESeq2(mat,cond)
  #dsres[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  tres <- ResDESeq2(mat,cond,sf=sfs)
  dsresq[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  #tres <- ResedgeR(mat,cond)
  #tres <-tres[rownames(mat),]
  #edres[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  tres <- ResedgeR(mat,cond,sf=sfs)
  tres <-tres[rownames(mat),]
  edresq[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  #tres <- ResVoom(mat,cond)
  #tres <-tres[rownames(mat),]
  #vmres[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  tres <- ResVoom(mat,cond,sf=sfs)
  tres <-tres[rownames(mat),]
  vmresq[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  #tres <- ResbaySeq(mat,cond)
  #bsres[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  #tres <- ResbaySeq(mat,cond,sf=sfs)
  #bsresq[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  tres <- ResABSSeq(mat,cond)
  abres[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  #not suit for ROTS
 # tres <- ResROTS(log2(totalnorm(mat)+1),cond)
 # rtres[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
 # tres <- ResROTS(log2(qtotalnorm(mat)+1),cond)
 # rtresq[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
}
basn="CT_5_"
path="./dataset/SimStudyCT_B_625_625/5samplespergroup"
cond<-c(rep("condA",5),rep("condB",5))
for(i in 1:10)
{
  filen<-paste(path,"/",basn,i,".rds",sep="")
  dat <-readRDS(filen)
  mat <-dat[[1]]
  gname <- rownames(mat)
  pset <- as.character(dat[[3]])
  indx <- gname%in%pset
  nset <- gname[!indx]
  beta <- rep(0,nrow(mat))
  beta[indx] <- log(1.5)
  sfs <- qtotalNormalized(mat)
  tres <- ResaFold(mat,cond)
  i <-i+10
  afres[i,] <-statsCal(tres[,1],adpv=tres[,3],beta=beta,gname=gname)
  tres <- ResaFold(mat,cond,normM = "TMM")
  #afrest[i,] <-statsCal(tres[,1],adpv=tres[,3],beta=beta,gname=gname)
  #tres <- ResaFold(mat,cond,normM = "geometric")
  #afresg[i,] <-statsCal(tres[,1],adpv=tres[,3],beta=beta,gname=gname)
  #tres <- ResaFold(mat,cond,normM = "quartile")
  #afresq[i,] <-statsCal(tres[,1],adpv=tres[,3],beta=beta,gname=gname)
  #tres <- ResDESeq2(mat,cond)
  #dsres[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  tres <- ResDESeq2(mat,cond,sf=sfs)
  dsresq[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
 # tres <- ResedgeR(mat,cond)
 # tres <-tres[rownames(mat),]
 # edres[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  tres <- ResedgeR(mat,cond,sf=sfs)
  tres <-tres[rownames(mat),]
  edresq[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  #tres <- ResVoom(mat,cond)
  #tres <-tres[rownames(mat),]
  #vmres[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  tres <- ResVoom(mat,cond,sf=sfs)
  tres <-tres[rownames(mat),]
  vmresq[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  #tres <- ResbaySeq(mat,cond)
  #bsres[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  #tres <- ResbaySeq(mat,cond,sf=sfs)
  #bsresq[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  tres <- ResABSSeq(mat,cond)
  abres[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  #tres <- ResROTS(log2(totalnorm(mat)+1),cond)
  #rtres[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  tres <- ResROTS(log2(qtotalnorm(mat)+1),cond)
  rtresq[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
}
basn="CT_10_"
path="./dataset/SimStudyCT_B_625_625/10samplespergroup"
cond<-c(rep("condA",10),rep("condB",10))
for(i in 1:10)
{
  filen<-paste(path,"/",basn,i,".rds",sep="")
  dat <-readRDS(filen)
  mat <-dat[[1]]
  gname <- rownames(mat)
  pset <- as.character(dat[[3]])
  indx <- gname%in%pset
  nset <- gname[!indx]
  beta <- rep(0,nrow(mat))
  beta[indx] <- log(1.5)
  sfs <- qtotalNormalized(mat)
  tres <- ResaFold(mat,cond)
  i <-i+20
  afres[i,] <-statsCal(tres[,1],adpv=tres[,3],beta=beta,gname=gname)
  #tres <- ResaFold(mat,cond,normM = "TMM")
  #afrest[i,] <-statsCal(tres[,1],adpv=tres[,3],beta=beta,gname=gname)
  #tres <- ResaFold(mat,cond,normM = "geometric")
  #afresg[i,] <-statsCal(tres[,1],adpv=tres[,3],beta=beta,gname=gname)
  #tres <- ResaFold(mat,cond,normM = "quartile")
  #afresq[i,] <-statsCal(tres[,1],adpv=tres[,3],beta=beta,gname=gname)
  #tres <- ResDESeq2(mat,cond)
  #dsres[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  tres <- ResDESeq2(mat,cond,sf=sfs)
  dsresq[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  #tres <- ResedgeR(mat,cond)
  #tres <-tres[rownames(mat),]
  #edres[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  tres <- ResedgeR(mat,cond,sf=sfs)
  tres <-tres[rownames(mat),]
  edresq[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  #tres <- ResVoom(mat,cond)
  #tres <-tres[rownames(mat),]
  #vmres[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  tres <- ResVoom(mat,cond,sf=sfs)
  tres <-tres[rownames(mat),]
  vmresq[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  #tres <- ResbaySeq(mat,cond)
  #bsres[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  #tres <- ResbaySeq(mat,cond,sf=sfs)
  #bsresq[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  tres <- ResABSSeq(mat,cond)
  abres[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  #tres <- ResROTS(log2(totalnorm(mat)+1),cond)
  #rtres[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
  tres <- ResROTS(log2(qtotalnorm(mat)+1),cond)
  rtresq[i,] <-statsCal(tres[,1],tres[,2],tres[,3],beta,gname)
}

NBres=list("aFold"=afres,"aFoldt"=afrest,"aFoldg"=afresg,"aFoldq"=afresq,
           "DESeq2"=dsres,"DESeq2q"=dsresq,"edgeR"=edres,"edgeRq"=edresq,
           "Voom"=vmres,"Voomq"=vmresq,"baySeq"=bsres,"baySeqq"=bsresq,
           "ABSSeq"=abres,"ROTS"=rtres,"ROTSq"=rtresq)
NBres[["aFold"]] <- afres
NBres[["DESeq2q"]] <- dsresq
NBres[["edgeRq"]] <- edresq
NBres[["Voomq"]] <- vmresq
NBres[["baySeqq"]] <- bsresq
NBres[["ABSSeq"]] <- abres
NBres[["ROTSq"]] <- rtresq
save(NBres,file="NB625res.Rdata")


load("NB625res.Rdata")
load("bsNB625res.Rdata")
NBres[["baySeq"]] <- bsNBres[["baySeq"]]
NBres[["baySeqq"]] <- bsNBres[["baySeqq"]]
save(NBres,file="NB625res.Rdata")

NBRres <- NBres
###boxplot
pdf("NB_625_625.pdf",width=12,height=8)

pd=cbind(NBRres[["aFold"]][1:10,3],NBRres[["DESeq2"]][1:10,3],NBRres[["edgeR"]][1:10,3],
         NBRres[["Voom"]][1:10,3],NBRres[["baySeq"]][1:10,3],NBRres[["ABSSeq"]][1:10,3],NBRres[["ROTS"]][1:10,3],
         NBRres[["aFold"]][11:20,3],NBRres[["DESeq2"]][11:20,3],NBRres[["edgeR"]][11:20,3],
         NBRres[["Voom"]][11:20,3],NBRres[["baySeq"]][11:20,3],NBRres[["ABSSeq"]][11:20,3],NBRres[["ROTS"]][11:20,3],
         NBRres[["aFold"]][21:30,3],NBRres[["DESeq2"]][21:30,3],NBRres[["edgeR"]][21:30,3],
         NBRres[["Voom"]][21:30,3],NBRres[["baySeq"]][21:30,3],NBRres[["ABSSeq"]][21:30,3],NBRres[["ROTS"]][21:30,3])
#colnames(pd)=rep(c("ABSSeq","edgeR","DESeq","baySeq","DESeq2","limmaVoom","limmaQN"),3)
boxplot(pd,at=c(1:7,9:15,17:23),frame=T,xlab="",ylab="",main="",axes=F,xlim=c(0.5,23.5),cex.axis=0.7,ylim=c(0.65,1),col=c(2:8))
axis(2,at=seq(0.7,1.0,0.1),labels=rep("",length(seq(0.7,1.0,0.1))))
#axis(1,at=1:7,labels=rep("",7))

abline(h=seq(0.7,1.0,.1),lty=2,col="gray")
abline(v=c(8,16),lty=2,col="gray")

dev.off()

pdf("NB_625_625_scat.pdf",width=12,height=8)

px=c(NBRres[["aFold"]][1:10,2],NBRres[["DESeq2"]][1:10,2],NBRres[["edgeR"]][1:10,2],
     NBRres[["Voom"]][1:10,2],NBRres[["baySeq"]][1:10,2],NBRres[["ABSSeq"]][1:10,2],NBRres[["ROTS"]][1:10,2],
     NBRres[["aFold"]][11:20,2]+1.1,NBRres[["DESeq2"]][11:20,2]+1.1,NBRres[["edgeR"]][11:20,2]+1.1,
     NBRres[["Voom"]][11:20,2]+1.1,NBRres[["baySeq"]][11:20,2]+1.1,NBRres[["ABSSeq"]][11:20,2]+1.1,NBRres[["ROTS"]][11:20,2]+1.1,
     NBRres[["aFold"]][21:30,2]+2.2,NBRres[["DESeq2"]][21:30,2]+2.2,NBRres[["edgeR"]][21:30,2]+2.2,
     NBRres[["Voom"]][21:30,2]+2.2,NBRres[["baySeq"]][21:30,2]+2.2,NBRres[["ABSSeq"]][21:30,2]+2.2,NBRres[["ROTS"]][21:30,2]+2.2)
py=c(NBRres[["aFold"]][1:10,1],NBRres[["DESeq2"]][1:10,1],NBRres[["edgeR"]][1:10,1],
     NBRres[["Voom"]][1:10,1],NBRres[["baySeq"]][1:10,1],NBRres[["ABSSeq"]][1:10,1],NBRres[["ROTS"]][1:10,2],
     NBRres[["aFold"]][11:20,1],NBRres[["DESeq2"]][11:20,1],NBRres[["edgeR"]][11:20,1],
     NBRres[["Voom"]][11:20,1],NBRres[["baySeq"]][11:20,1],NBRres[["ABSSeq"]][11:20,1],NBRres[["ROTS"]][11:20,1],
     NBRres[["aFold"]][21:30,1],NBRres[["DESeq2"]][21:30,1],NBRres[["edgeR"]][21:30,1],
     NBRres[["Voom"]][21:30,1],NBRres[["baySeq"]][21:30,1],NBRres[["ABSSeq"]][21:30,1],NBRres[["ROTS"]][21:30,1])
#colnames(pd)=rep(c("ABSSeq","edgeR","DESeq","baySeq","DESeq2","limmaVoom","limmaQN"),3)
plot(px,py,frame=T,xlab="",ylab="",main="",axes=F,xlim=c(0,3.2),cex.axis=0.7,ylim=c(0,1),col=rep(c(2:8),each=10),pch=16)
axis(2,at=seq(0.0,1.0,0.25),labels=F)
axis(1,at=c(seq(0,1,0.25),seq(1.1,2.1,0.25),seq(2.2,3.2,0.25)),labels=F)
#rect(12,0.840,13,0.848,col=2)
#rect(12,0.819,13,0.827,col=3)
#rect(12,0.798,13,0.806,col=4)
#rect(12,0.777,13,0.785,col=5)
#rect(12,0.756,13,0.764,col=6)
#rect(12,0.735,13,0.743,col=7)

abline(h=seq(0.,1.0,.25),lty=2,col="gray")
abline(v=c(1.05,2.15),lty=1,col="black")
abline(v=c(.1,1.2,2.3),lty=2,col="gray")
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



pdf("NB-q_625_625.pdf",width=12,height=8)

pd=cbind(NBRres[["aFold"]][1:10,3],NBRres[["aFoldt"]][1:10,3],NBRres[["aFoldg"]][1:10,3],NBRres[["aFoldq"]][1:10,3]
         ,NBRres[["DESeq2"]][1:10,3],NBRres[["DESeq2q"]][1:10,3],NBRres[["edgeR"]][1:10,3],NBRres[["edgeRq"]][1:10,3],
         NBRres[["Voom"]][1:10,3],NBRres[["Voomq"]][1:10,3],NBRres[["baySeq"]][1:10,3],NBRres[["ABSSeq"]][1:10,3],NBRres[["ROTS"]][1:10,3],
         NBRres[["aFold"]][11:20,3],NBRres[["aFoldt"]][11:20,3],NBRres[["aFoldg"]][11:20,3],NBRres[["aFoldq"]][11:20,3],
         NBRres[["DESeq2"]][11:20,3],NBRres[["edgeR"]][11:20,3],
         NBRres[["Voom"]][11:20,3],NBRres[["baySeq"]][11:20,3],NBRres[["ABSSeq"]][11:20,3],NBRres[["ROTS"]][11:20,3],
         NBRres[["aFold"]][21:30,3],NBRres[["aFoldt"]][21:30,3],NBRres[["aFoldg"]][21:30,3],NBRres[["aFoldq"]][21:30,3],
         NBRres[["DESeq2"]][21:30,3],NBRres[["edgeR"]][21:30,3],
         NBRres[["Voom"]][21:30,3],NBRres[["baySeq"]][21:30,3],NBRres[["ABSSeq"]][21:30,3],NBRres[["ROTS"]][21:30,3])
#colnames(pd)=rep(c("ABSSeq","edgeR","DESeq","baySeq","DESeq2","limmaVoom","limmaQN"),3)
boxplot(pd,at=c(1:7,9:15,17:23),frame=T,xlab="",ylab="",main="",axes=F,xlim=c(0.5,23.5),cex.axis=0.7,ylim=c(0.65,1),col=c(2:8))
axis(2,at=seq(0.7,1.0,0.1),labels=rep("",length(seq(0.7,1.0,0.1))))
#axis(1,at=1:7,labels=rep("",7))

abline(h=seq(0.7,1.0,.1),lty=2,col="gray")
abline(v=c(8,16),lty=2,col="gray")

dev.off()

pdf("NB-q_625_625_scat.pdf",width=12,height=8)

px=c(NBRres[["aFold"]][1:10,2],NBRres[["DESeq2"]][1:10,2],NBRres[["edgeR"]][1:10,2],
     NBRres[["Voom"]][1:10,2],NBRres[["baySeq"]][1:10,2],NBRres[["ABSSeq"]][1:10,2],NBRres[["ROTS"]][1:10,1],
     NBRres[["aFold"]][11:20,2]+1.1,NBRres[["DESeq2"]][11:20,2]+1.1,NBRres[["edgeR"]][11:20,2]+1.1,
     NBRres[["Voom"]][11:20,2]+1.1,NBRres[["baySeq"]][11:20,2]+1.1,NBRres[["ABSSeq"]][11:20,2]+1.1,NBRres[["ROTS"]][11:20,2]+1.1,
     NBRres[["aFold"]][21:30,2]+2.2,NBRres[["DESeq2"]][21:30,2]+2.2,NBRres[["edgeR"]][21:30,2]+2.2,
     NBRres[["Voom"]][21:30,2]+2.2,NBRres[["baySeq"]][21:30,2]+2.2,NBRres[["ABSSeq"]][21:30,2]+2.2,NBRres[["ROTS"]][21:30,2]+2.2)
py=c(NBRres[["aFold"]][1:10,1],NBRres[["DESeq2"]][1:10,1],NBRres[["edgeR"]][1:10,1],
     NBRres[["Voom"]][1:10,1],NBRres[["baySeq"]][1:10,1],NBRres[["ABSSeq"]][1:10,1],NBRres[["ROTS"]][1:10,1],
     NBRres[["aFold"]][11:20,1],NBRres[["DESeq2"]][11:20,1],NBRres[["edgeR"]][11:20,1],
     NBRres[["Voom"]][11:20,1],NBRres[["baySeq"]][11:20,1],NBRres[["ABSSeq"]][11:20,1],NBRres[["ROTS"]][11:20,1],
     NBRres[["aFold"]][21:30,1],NBRres[["DESeq2"]][21:30,1],NBRres[["edgeR"]][21:30,1],
     NBRres[["Voom"]][21:30,1],NBRres[["baySeq"]][21:30,1],NBRres[["ABSSeq"]][21:30,1],NBRres[["ROTS"]][21:30,1])

apx=c(NBRres[["aFoldt"]][1:10,2],NBRres[["aFoldg"]][1:10,2],NBRres[["aFoldq"]][1:10,2],
     NBRres[["aFoldt"]][11:20,2]+1.1,NBRres[["aFoldg"]][11:20,2]+1.1,NBRres[["aFoldq"]][11:20,2]+1.1,
     NBRres[["aFoldt"]][21:30,2]+2.2,NBRres[["aFoldg"]][21:30,2]+2.2,NBRres[["aFoldq"]][21:30,2]+2.2)
apy=c(NBRres[["aFoldt"]][1:10,1],NBRres[["aFoldg"]][1:10,1],NBRres[["aFoldq"]][1:10,1],
     NBRres[["aFoldt"]][11:20,1],NBRres[["aFoldg"]][11:20,1],NBRres[["aFoldq"]][11:20,1],
     NBRres[["aFoldt"]][21:30,1],NBRres[["aFoldg"]][21:30,1],NBRres[["aFoldq"]][21:30,1])
bpx=c(NBRres[["DESeq2q"]][1:10,2],NBRres[["edgeRq"]][1:10,2],
     NBRres[["Voomq"]][1:10,2],NBRres[["baySeqq"]][1:10,2],NBRres[["ROTSq"]][1:10,1],
     NBRres[["DESeq2q"]][11:20,2]+1.1,NBRres[["edgeRq"]][11:20,2]+1.1,
     NBRres[["Voomq"]][11:20,2]+1.1,NBRres[["baySeqq"]][11:20,2]+1.1,NBRres[["ROTSq"]][11:20,2]+1.1,
     NBRres[["DESeq2q"]][21:30,2]+2.2,NBRres[["edgeRq"]][21:30,2]+2.2,
     NBRres[["Voomq"]][21:30,2]+2.2,NBRres[["baySeqq"]][21:30,2]+2.2,NBRres[["ROTSq"]][21:30,2]+2.2)
bpy=c(NBRres[["DESeq2q"]][1:10,1],NBRres[["edgeRq"]][1:10,1],
     NBRres[["Voomq"]][1:10,1],NBRres[["baySeqq"]][1:10,1],NBRres[["ROTSq"]][1:10,1],
     NBRres[["DESeq2q"]][11:20,1],NBRres[["edgeRq"]][11:20,1],
     NBRres[["Voomq"]][11:20,1],NBRres[["baySeqq"]][11:20,1],NBRres[["ROTSq"]][11:20,1],
     NBRres[["DESeq2q"]][21:30,1],NBRres[["edgeRq"]][21:30,1],
     NBRres[["Voomq"]][21:30,1],NBRres[["baySeqq"]][21:30,1],NBRres[["ROTSq"]][21:30,1])

#colnames(pd)=rep(c("ABSSeq","edgeR","DESeq","baySeq","DESeq2","limmaVoom","limmaQN"),3)
plot(px,py,frame=T,xlab="",ylab="",main="",axes=F,xlim=c(0,3.2),cex.axis=0.7,ylim=c(0,1),col=rep(c(2:8),each=10),pch=16)
points(apx,apy,col=2,pch=c(3:5))
points(bpx,bpy,col=rep(c(3:6,8),each=10),pch=1)
axis(2,at=seq(0.0,1.0,0.25),labels=F)
axis(1,at=c(seq(0,1,0.25),seq(1.1,2.1,0.25),seq(2.2,3.2,0.25)),labels=F)

points(0.12,0.95,pch=1)
points(0.12,0.85,pch=3)
points(0.12,0.75,pch=4)
points(0.12,0.65,pch=5)

abline(h=seq(0.,1.0,.25),lty=2,col="gray")
abline(v=c(1.05,2.15),lty=1,col="black")
abline(v=c(.1,1.2,2.3),lty=2,col="gray")
dev.off()