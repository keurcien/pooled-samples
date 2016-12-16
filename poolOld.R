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
  res$pvalues <- compute.pval(res$chi2.stat,K,method="mahalanobis")
  class(res) <- 'pcadapt'
  attr(res,"K") <- K
  return(res)
}

pool.old.corrected = function(data,K,min.maf,cover.matrix=NULL){
  nSNP <- ncol(data)
  nPOP <- nrow(data)
  # New procedure
  z.matrix <- array(0,dim=c(nPOP,nSNP))
  for (k in 1:nSNP){
    for (n in 1:nPOP){
      n_i <- cover.matrix[n,k]
      f_i <- data[n,k]
      se <- sqrt(abs(f_i*(1-f_i))/n_i)
      if ((!is.na(se)) && (se > 0)){
        z.matrix[n,k] <- f_i/se 
      } else {
        z.matrix[n,k] <- f_i
      }
    }
  }
  # End new procedure
  res <- corpca(data=data,K=K)
  freq <- apply(data,2,FUN=function(x){mean(x,na.rm=TRUE)})
  res$maf <- as.vector(pmin(freq,1-freq))
  res$loadings[res$maf<min.maf] <- NA 
  z.matrix[,res$maf<min.maf] <- NA
  res$stat <- array(NA,dim=nSNP)
  finite.list <- which(!is.na(apply(abs(z.matrix),2,sum)))
  res$stat <- as.vector(robust::covRob(t(z.matrix),na.action=na.omit,estim = "pairwiseGK")$dist)
  res$chi2.stat <- res$stat
  # Compute p-values
  res$pvalues <- compute.pval(res$chi2.stat,K,method="mahalanobis")
  class(res) <- 'pcadapt'
  attr(res,"K") <- K
  attr(res,"method") <- "mahalanobis"
  attr(res,"data.type") <- "pool"
  attr(res,"min.maf") <- min.maf
  return(res)
}