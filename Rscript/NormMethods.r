library(ABSSeq)
library(cqn)
library(edgeR)
library(DESeq2)


##Med-pgQ2 and UQ-pgQ2; from Li et. al. (2017) http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0176185

uq<-function(X){
  
  #excluding zero counts in each sample
  UQ<-function(y){
    quantile(y, 0.75)
  }
  X<-X+0.1
  upperQ<-apply(X,2,UQ)
  f.uq<-upperQ/mean(upperQ)
  upq.res<-scale(X,center=FALSE,scale=f.uq) 
  return(upq.res)
}

## Med normalization
med<-function(X){
  
  MED<-function(y){
    median(y[y>0])
  }
  X<-X+0.1
  med<-apply(X,2, MED)
  f.med<-med/mean(med)
  med.res<-scale(X,center=FALSE, scale=f.med)
  return(med.res)
}

# per gene normalization by Median: pgQ2
# X: a matrix of data with the multiplication of factor (f) as:
# f=50, 100, 200, 500 or 1000

# run with f=100 as suggested by Li et al. (2017)
pgene1<-function(X, f=100){
  m<-apply(X,1, median)
  si<-m/f # multiply f=100 per gene per sample
  X1<-scale(t(X),center=FALSE,scale=si)
  res<-t(X1)
  rownames(res)<-rownames(X)
  return(res)  
}
#MedpgQ2
MedpgQ2=function(x,f=100)
{
  norm.med<-med(x)
  return(pgene1(norm.med, f))
}
#UQpgQ2
UQpgQ2=function(x,f=100)
{
  norm.med<-uq(x)
  return(pgene1(norm.med, f))
}
#CQN
cqnnorm=function(x,gc,len)
{
  buf <-cqn(x,gc,len)
  return(buf$y+buf$offset)
}
#qtotal
qtotalnorm1=function(x)
{
  sf<-ABSSeq::qtotalNormalized(x)
  sf=sf/mean(sf)
  return( t( t(x) / sf) )
}
#qtotal
qtotalnorm=function(x)
{
  sf<-qtotalNormalized(x)
  sf=sf/mean(sf)
  return( t( t(x) / sf) )
}
#total
totalnorm=function(x)
{
  sf<-colSums(x)
  sf<-sf/mean(sf)
  return( t( t(x) / sf) )
}
#QN baySeq
upperqnorm=function(x)
{
  rowQuar=function(z) {
    x <- z[z > 0]
    sum(x[x <= quantile(x, 0.75, na.rm = TRUE)], na.rm = TRUE)
  }
  cbuf <- x
  cbuf <- cbuf[colSums(cbuf)>0,]
  sf <- apply(cbuf,2,rowQuar)
  sf<-sf/mean(sf)
  return( t( t(x) / sf) )
}

#DESEQ & DESeq2
deseqnorm=function(x)
{
  sf <- estimateSizeFactorsForMatrix(x)
  sf<-sf/mean(sf)
  return( t( t(x) / sf) )
}

#TMM
tmmnorm=function(x)
{
  nf <- calcNormFactors(x)
  sf<-colSums(x) * nf
  sf<-sf/mean(sf)
  return( t( t(x) / sf) )
}
