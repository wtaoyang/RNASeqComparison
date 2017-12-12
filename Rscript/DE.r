##aFold
library(ABSSeq)
ResaFold = function(x,cond,normM = "qtotal",sf)
{
  cds = ABSDataSet(x, cond,normMethod = normM)
  if(!missing(sf))  cds = ABSDataSet(x, cond,normMethod = "user",sizeFactor = sf)
  res=ABSSeq(cds,useaFold =T)
  res5=ABSSeq::results(res,c("foldChange","pvalue","adj.pvalue"))
  return(res5)
}

##ABSSeq
ResABSSeq = function(x,cond,normM = "qtotal")
{
  cds = ABSDataSet(x, cond,normMethod = normM)
  res=ABSSeq(cds,useaFold =F)
  res5=ABSSeq::results(res,c("foldChange","pvalue","adj.pvalue"))
  return(res5)
}
library(DESeq2)
#sf sizefactor
ResDESeq2=function(x,cond,sf)
{
  
  DESeq2.ds = DESeq2::DESeqDataSetFromMatrix(countData = x, colData = data.frame(cond), design = ~ cond)
  if(!missing(sf)) sizeFactors(DESeq2.ds)<-sf
  DESeq2.ds = DESeq2::DESeq(DESeq2.ds)
  
  DESeq2.results <- DESeq2::results(DESeq2.ds)
  DESeq2.results <- lfcShrink(dds=DESeq2.ds, contrast=c("cond","condB","condA"), res=DESeq2.results)
  DESeq2.pvalues <- DESeq2.results$pvalue 
  DESeq2.adjpvalues <- DESeq2.results$padj
  DESeq2.logFC <- DESeq2.results$log2FoldChange
  result.table <- cbind(DESeq2.logFC,DESeq2.pvalues,DESeq2.adjpvalues)
  rownames(result.table) <- rownames(x)
  colnames(result.table) <- c("logFC","pvalue","advalue")
  naindx <-is.na(result.table[,"advalue"])
  result.table[naindx,"advalue"] <- 1
  result.table[naindx,"pvalue"] <- 1
  return(result.table)
}

#sf sizefactor
library(limma)
library(edgeR)
ResVoom <- function(x, cond,sf,q.cut=0.05, lfc=0.0, condA='condA', condB='condB'){
  groups=as.factor(cond)
  design <- model.matrix(~0+groups)
  colnames(design) <- levels(groups)
  count.dat<-x
  lib.size <- 1
  if(!missing(sf)) lib.size <- sf
  else
  {
    nf <- calcNormFactors(count.dat)
    lib.size<-colSums(count.dat) * nf
  }
  dat <- voom(count.dat, design, plot=FALSE, lib.size=lib.size)
  fit=lmFit(dat,design)
  contrast.matrix <- makeContrasts("condB - condA", levels=design)
  
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2)
  res=decideTests(fit2,p.value=q.cut,lfc=lfc)
  tab<-topTable(fit2, adjust = "BH", number=nrow(fit2), sort.by='logFC')
  tab <- tab[,c("logFC","P.Value","adj.P.Val")]
  return(tab)
}

##edgeR
#sf sizefactor
library(edgeR)
ResedgeR <- function(x, cond,sf){
  y <- c()
  if(!missing(sf))
  {
    y <- DGEList(counts=x, group=factor(cond),lib.size=sf)
    y <- calcNormFactors(y,method="none")
  }else
  {
    y <- DGEList(counts=x, group=factor(cond))
    y <- calcNormFactors(y)
  }
  ## estimate common dispersion
  y <- estimateCommonDisp(y)
  ## estimate gene specific dispersion
  y <- estimateTagwiseDisp(y,trend="movingave")
  ## DE test
  et <- exactTest(y, dispersion="tagwise")
  res <- topTags(et,n=nrow(et))[[1]]
  res <- res[,c(1,3,4)]
  return(res)
}


library("baySeq")
library("snow")
#cl number of Threads
ResbaySeq=function(x,cond,sf,cl=NULL)
{
  if(!is.null(cl))
  {
    cl <-makeCluster(cl, "SOCK")
  }
  baySeq.cd <- new("countData", data = as.matrix(x), replicates = cond, groups = list(NDE = rep(1, length(cond)), DE = cond))
  libsizes(baySeq.cd) <- baySeq::getLibsizes(baySeq.cd)
  if(!missing(sf)) libsizes(baySeq.cd) <- sf
  baySeq.cd <- baySeq::getPriors.NB(baySeq.cd,smaplesize=5000, cl = cl) 
  baySeq.cd <- baySeq::getLikelihoods(baySeq.cd, cl = cl)
  stopCluster(cl)
  baySeq.posteriors.DE <- exp(baySeq.cd@posteriors)[, 2] 
  baySeq.FDR <- baySeq::topCounts(baySeq.cd, group = 'DE', FDR = 1)$FDR[match(rownames(dat[[1]]), rownames(baySeq::topCounts(baySeq.cd, group = 'DE', FDR = 1)))] 
  result.table <- cbind(baySeq.posteriors.DE,baySeq.FDR) 
  rownames(result.table) = rownames(x)
  colnames(result.table) = c("DEsign","FDR")
  return(result.table)
}

library("ROTS")
#x mast be positive in order to log transformed
ResROTS=function(x,cond)
{
  results = ROTS(data = x, groups = cond , B = 100 , K = 500 , seed = 1234)
  result.table <- cbind(results$logfc,results$pvalue,results$FDR)
  return(result.table)
}
