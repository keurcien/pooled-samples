#' Convert genotypes to pooled samples
#'
#' \code{get.pool.matrix} creates a pooled-sequenced data out of a genotype matrix,
#' given the labels of each individuals.
#' 
#' @param data a matrix with n rows and p columns where n is the number of individuals and p is the number of markers. 
#' @param pop a list of integers or strings specifying which subpopulation the individuals belong to.
#' @param ploidy an integer specifying the ploidy of the individuals.
#'
get.pool.matrix = function(data,pop,ploidy=2){
  nSNP <- ncol(data)
  pop.names <- unique(pop)
  freq <- array(0,dim=c(length(pop.names),nSNP))
  for (k in 1:length(pop.names)){
    geno.k <- data[pop==pop.names[k],]
    geno.k[geno.k==9] <- NA
    freq[k,] <- apply(geno.k,MARGIN=2,FUN=function(xx){mean(xx,na.rm=TRUE)})/ploidy
  }
  return(freq)
}

#' Simulate coverage matrix
#'
#' \code{simulate.cover} creates a matrix simulating the coverage by SNP and by pool. 
#' 
#' @param nrow an integer equal to the number of pools.
#' @param ncol an integer equal to the number of markers.
#' @param min.cover an integer indicating the minimum coverage.
#' @param max.cover an integer indicating the maximum coverage.
#' @param high.cov.loc a list of integers indicating the indices of SNPs with higher coverage.
#' @param high.cov.pool a list of integers indicating the indices of pools with globally higher coverage.
#'
simulate.cover = function(nrow,ncol,min.cover,max.cover,high.cov.loc=NULL,high.cov.pool=NULL){
  if (length(min.cover)==1){
    rnum <- floor(runif(n = nrow*ncol,min = min.cover,max = max.cover))+1
    coverage.mat <- matrix(rnum,nrow = nrow,ncol = ncol)
  } else if (length(min.cover)==2){
    if (!is.null(high.cov.loc)){
      hcl <- high.cov.loc
    } else {
      if (ncol>20){
        hcl <- (ncol-20):(ncol)
      } else {
        hcl <- (ncol-1):(ncol)
      }
    }
    if (is.null(high.cov.pool)){
      n.hc <- length(hcl)
      rnum <- floor(runif(n = nrow*ncol,min = min.cover[1],max = max.cover[1]))+1  
      rnum.hc <- floor(runif(n = nrow*n.hc,min = min.cover[2],max = max.cover[2]))+1
      hc.mat <- matrix(rnum.hc,nrow = nrow,ncol = n.hc)
      coverage.mat <- matrix(rnum,nrow = nrow,ncol = ncol)
      coverage.mat[,hcl] <- hc.mat
    } else {
      n.hcp <- length(high.cov.pool)
      n.hc <- length(hcl)
      rnum <- floor(runif(n = nrow*ncol,min = min.cover[1],max = max.cover[1]))+1  
      rnum.hc <- floor(runif(n = n.hcp*n.hc,min = min.cover[2],max = max.cover[2]))+1
      hc.mat <- matrix(rnum.hc,nrow = n.hcp,ncol = n.hc)
      coverage.mat <- matrix(rnum,nrow = nrow,ncol = ncol)
      coverage.mat[high.cov.pool,hcl] <- hc.mat  
    }
  } else {
    stop("Unvalid min.cover.")
  }
  return(coverage.mat)
}


#' Sample genotype matrix from pooled samples
#'
#' \code{sample.geno} samples a genotype matrix from pooled samples.
#' 
#' @param pool.matrix a matrix with n rows and p columns where n is the number of pools and is the number of markers.
#' @param cover.matrix a matrix with n rows and p columns where n is the number of pools and is the number of markers.
#' @param nINDperPOP a list specifying the number of individuals for each pool.
#' @param method a character string indicating the method used for sampling.
#'
sample.geno = function(pool.matrix=NULL,cover.matrix=NULL,nINDperPOOL=NULL,method="betabino"){
  nPOOL <- nrow(pool.matrix)  
  nSNP <- ncol(pool.matrix)
  if (missing(nINDperPOOL)){
    sample.size <- rep(20,nPOOL)
  } else {
    sample.size <- nINDperPOOL
  }
  if (missing(cover.matrix)){
    print("Coverage matrix missing. Drawing genotypes from a binomial distribution.")
    geno <- NULL
    for (k in 1:nPOOL){
      G <- rbinom(n=(sample.size[k]*nSNP),size=2,prob=pool.matrix[k,])
      geno <- rbind(geno,t(matrix(G,nrow = nSNP,sample.size[k])))
    }
  } else {
    print("Coverage matrix not missing. Drawing genotypes from a beta-binomial distribution.")
    geno <- NULL
    for (k in 1:nPOOL){
      cover.1 <- cover.matrix[k,]
      n.reads <- array(NA,dim=c(1,nSNP))
      nna <- which(pool.matrix[k,]>0)
      na <- which(pool.matrix[k,]==0)
      epsilon <- 0.00001
      n.reads[nna] <- cover.1[nna]/pool.matrix[k,nna]
      n.reads[na] <- cover.1[na]/epsilon
      cover.2 <- n.reads - cover.1
      if (method=="betabino"){
        p <- matrix(rbeta(sample.size[k]*nSNP,cover.1+1,cover.2+1),nrow=sample.size[k],ncol=nSNP)
        p.aux <- matrix(p,nrow=1)
        G <- t(matrix(rbinom(sample.size[k]*nSNP,size=2,prob = p.aux),ncol=sample.size[k],nrow=nSNP))
        geno <- rbind(geno,G)
      } 
#      geno <- rbind(geno,G)
#       for (i in 1:sample.size[k]){
#         p <- array(0,dim=c(1,nSNP))
#         p[,] <- sapply(1:nSNP,FUN=function(l){rbeta(n = 1,cover.1[l]+1,cover.2[l]+1)})
#         G <- array(0,dim=c(1,nSNP))
#         G[,] <- sapply(1:nSNP,FUN=function(l){rbinom(n = 1, size = 2, prob = p[1,l])})
#         geno <- rbind(geno,G)
#       }
    }
  }
  return(geno)
}

cover.to.pool = function(data,cover.matrix,pop,ploidy=2){
  nPOP <- nrow(cover.matrix)
  nSNP <- ncol(cover.matrix)
  pool.matrix <- array(0,dim=c(nPOP,nSNP))
  pop.lab <- unique(pop)
  for (n in 1:nPOP){
    idx <- which(pop==pop.lab[n])
    nIND.pop.n <- length(idx)
    for (p in 1:nSNP){
      c.np <- cover.matrix[n,p]
      if (c.np > 0){
        draw.np <- floor(runif(c.np,min = 1,max = nIND.pop.n) + 1)
        drawn.geno <- data[idx[draw.np],p]
        pool.matrix[n,p] <- sum(drawn.geno,na.rm = TRUE)/(ploidy*c.np)
      } else {
        pool.matrix[n,p] <- NA
      }
    }
  }
  return(pool.matrix)
}

