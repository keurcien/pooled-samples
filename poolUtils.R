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
#'
simulate.cover = function(nrow,ncol,min.cover,max.cover){
  rnum <- floor(runif(n = nrow*ncol,min=min.cover,max=max.cover))+1
  return(matrix(rnum,nrow = nrow,ncol = ncol))
}


#' Sample genotype matrix from pooled samples
#'
#' \code{sample.geno} samples a genotype matrix from pooled samples.
#' 
#' @param pool.matrix a matrix with n rows and p columns where n is the number of pools and is the number of markers.
#' @param cover.matrix a matrix with n rows and p columns where n is the number of pools and is the number of markers.
#' @param nINDperPOP a list specifying the number of individuals for each pool.
#'
sample.geno = function(pool.matrix,cover.matrix,nINDperPOOL=NULL){
  nPOOL <- nrow(pool.matrix)  
  nSNP <- ncol(pool.matrix)
  if (missing(nINDperPOOL)){
    sample.size <- rep(nPOOL,20)
  } else {
    sample.size <- nINDperPOOL
  }
  geno <- NULL
  for (k in 1:nPOOL){
    one <- matrix(1,nrow= sample.size,1)
    rep.G[] <- one %*% pool.matrix[k,]
    G <- rep.G %>%
      map_dbl( ~ rbinom(n = 1, size = 2, prob = .x)) %>%
      matrix(nrow = d * pop.nbindiv, ncol = L)
    geno <- rbind(geno,G)
  }

  return(G)
}


