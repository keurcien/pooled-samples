library(pcadapt)
library(robust)
library(MASS)

pool.old <- function(input,K,min.maf=0.05){
  nPOP <- nrow(input)
  nSNP <- ncol(input)
  if (missing(K)){
    K <- nPOP-1
  }
  res <- corpca(input,K)
  freq <- apply(input,2,FUN=function(x){mean(x,na.rm=TRUE)})
  res$maf <- as.vector(pmin(freq,1-freq))
  res$loadings[res$maf<min.maf] <- NA 
  res$stat <- array(NA,dim=nSNP)
  finite.list <- which(!is.na(apply(abs(res$loadings),1,sum)))
  if (K>1){
    res$stat[finite.list] <- as.vector(robust::covRob(res$loadings,na.action=na.omit,estim="pairwiseGK")$dist)
  } else {
    onedcov <- as.vector(MASS::cov.rob(res$loadings[finite.list,1]))
    res$stat <- (res$zscores[,1]-onedcov$center)^2/onedcov$cov[1]
  }
  res$gif <- median(res$stat,na.rm=TRUE)/qchisq(0.5,df=K)
  res$chi2.stat <- res$stat/res$gif
  return(res)
}

pool.old.corrected <- function(input,K,min.maf=0.05){
  
}