##AUC calculator
library(ROC)
area<-function(x,y)
{
  x <- 1 - x
  if (x[1] > x[length(x)]) {
    x <- rev(x)
    y <- rev(y)
  }
  trapezint(x, y, 0, 1)
} 


#gname: ordered gene name
#pset: positive set, gene name
#nset: negative set, gene name
#return value: a list with x, y for curve and AUC value
roccurve=function(gname,pset,nset)
{
  x <- c()
  y <- c()
  z <- 0
  p <- length(pset)
  n <-  length(nset)
  for(i in seq(50,length(gname),50))
  {
    tp<-sum(pset%in%gname[1:i])
    fp<-sum(nset%in%gname[1:i])
    y <- c(y,tp/p)
    x <- c(x,fp/n)
  }
  z <- area(x,y)
  return(list(x=x,y=y,auc=z))
}

#gname: ordered gene name
#pset: positive set, gene name
#nset: negative set, gene name
#return value: a list with x, y for fdr sensitive curve
fdrcurve=function(gname,pset,nset)
{
  x <- c()
  y <- c()
  z <- 0
  p <- length(pset)
  n <-  length(nset)
  for(i in seq(50,length(gname),50))
  {
    tp<-sum(pset%in%gname[1:i])
    fp<-sum(nset%in%gname[1:i])
    y <- c(y,fp/(tp+fp))
    x <- c(x,i)
  }
  return(list(x=x,y=y))
}

#gname: ordered gene name
#pset: positive set, gene name
#nset: negative set, gene name
#return value: a list with x, y for fdr sensitive curve
fncurve=function(gname,pset,nset)
{
  x <- c()
  y <- c()
  z <- 0
  p <- length(pset)
  n <-  length(nset)
  for(i in seq(50,length(gname),50))
  {
    tp<-sum(pset%in%gname[1:i])
    fp<-sum(nset%in%gname[1:i])
    y <- c(y,fp)
    x <- c(x,i)
  }
  return(list(x=x,y=y))
}
#### A function for calculating the  standard error (se) of AUC values from Li et. al. (2017) 
#### http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0176185
# A: AUC value
# na: number of true positive genes
# nn: number of  true negative genes
se<-function(A, na, nn){
  Q1<-A/(2-A)
  Q2<-2*A^2/(1+A)
  d1<-A*(1-A)+(na-1)*(Q1-A^2)+(nn-1)*(Q2-A^2)
  S<-sqrt(d1/(na*nn))
  return(S)  
}