# Get frankensigs takes an input of a matrix of signatures, where the rows are features and the
# columns are signatures.  It generates frankensigs by averaging n signatures together, where
# n is the given parameter (default 1). 
# If return2 == 1, it returns a list of two disjoint franken signatures of size n. 
# If returnk == 1, it returns the maximum number of disjoint franken signatures of size n floor(ncols / n).
# If return2 == 0, it returns a single signature. 

get_frankensigs <- function(mymat, k=1, return2=1, returnk=0){
  if (k==1){
    if (return2 == 0){
      return(mymat[, sample(dim(mymat)[2], 1)])
    } else {
      ix <- sample(dim(mymat)[2], 2)
      return(cbind(mymat[, ix[1]], mymat[, ix[2]]))
    }
  }
  
  if (return2 == 0){
    ix <- sample(dim(mymat)[2], k)
    return(rowMeans(mymat[, ix]))
  } else if (returnk==1) {
    
    ix <- sample(dim(mymat)[2], floor(dim(mymat)[2]/k) * k)
    return(sapply(seq(floor(dim(mymat)[2]/k)), FUN=function(x) rowMeans(mymat[, ix[(k*(x-1)+1):(k*x)]])))

  } else {
    if (2*k <= dim(mymat)[2]){
      ix <- sample(dim(mymat)[2], 2*k)
      return(sapply(seq(2), FUN=function(x) rowMeans(mymat[, ix[(k*(x-1)+1):(k*x)]])))
      #return(list(rowMeans(mymat[, ix[1:n]]), rowMeans(mymat[, ix[(n+1):(2*n)]])))
    } else {
      return(c(0,0))
    }
  }
  
}



get_frankencorr <- function(mymat, nmax, iter=100, useall=FALSE){
  nmax <- min(nmax, floor(dim(mymat)[2]/2))
  
  # Revisit:
  # if (nmax > 50){
  #   seqvals <- c(seq(10), seq(15, nmax, 5))
  #   seqvals <- seqvals[seqvals <= nmax]
  # } else {
  #   seqvals <- seq(nmax)
  # }
  
  # if (!useall){
  #   seqvals <- c(seq(10), seq(15, nmax, 5))
  #   seqvals <- seqvals[seqvals <= nmax]
  # } else {
  #   seqvals <- seq(1, nmax)
  # }
  seqvals <- seq(100, nmax, 100)
  
  cormat <- matrix(numeric(iter*length(seqvals)), nrow=length(seqvals), dimnames=list(seqvals))
  
  for (jj in seq(length(seqvals))){
    ii <- seqvals[jj]
    a <- lapply(seq(iter), FUN=function(x) get_frankensigs(mymat, k=ii, return2=1))
    cormat[jj,] <- sapply(a, FUN=function(x) cor(x)[2])
    #cormat[ii,] <- sapply(seq(iter), FUN=function(x) cor(a[[1,x]], a[[2,x]]))
  }
  
  cordf <- data.frame(frankensize=seqvals, meancor=rowMeans(cormat), sdcor=apply(cormat,1,sd))
  return(list(cormat=cormat, cordf=cordf))
}


bootstrap_maxfrankencorr <- function(mymat, iter=100){
  nmax <- floor(dim(mymat)[2]/2)
  
  mycors <- sapply(seq(iter), FUN=function(x) cor(get_frankensigs(mymat, k=nmax, return2=1))[2])
  return(mycors)
}


# This function is inordinately slow
# This function takes as input a matrix and computes the 
get_frankensim <- function(mymat, nmax, iter=100, metric="pearson"){
  nmax <- min(nmax, floor(dim(mymat)[2]/2))
  
  if (nmax > 50){
    seqvals <- seq(5, nmax, 5)
  } else {
    seqvals <- seq(nmax)
  }
  
  simmat <- matrix(numeric(length(seqvals)*iter), nrow=length(seqvals))
  
  myfunc <- switch(metric, 
                   pearson = function(x) cor(x)[2], 
                   spearman = function(x) cor(x, method="spearman")[2], 
                   wtcs = function(x) compute_sim_block(GCT(mat=x, cid=as.character(seq(dim(x)[2]))), GCT(mat=x, cid=as.character(seq(dim(x)[2]))), metric=metric)[2], 
                   cosine = function(x) cosine(x)[2], 
                   wtcsmat = function(x) compute_cs_mats(x, x)[2])
  
  for (ii in seq_along(seqvals)){
    a <- lapply(seq(iter), FUN=function(x) get_frankensigs(mymat, k=seqvals[ii], return2=1))
    simmat[ii,] <- sapply(a, FUN=myfunc)
  }
  
  cordf <- data.frame(frankensize=seqvals, meancor=rowMeans(simmat), sdcor=apply(simmat,1,sd))
  rownames(simmat) <- seqvals
  return(list(simmat=simmat, cordf=cordf, metric=metric))
}


get_metasigs <- function(dspath, rid, siginfo, pert_iname=c()){
  retds <- c()
  for (mypert in pert_iname){
    mysigs <- siginfo$sig_id[siginfo$pert_iname == mypert]
    ds <- parse_gctx(dspath, rid=geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], cid=mysigs)
  }
}


ecdf_pointwise <- function(mydat, breaks){
  h <- hist(mydat, breaks=breaks, plot=FALSE)
  y0 <- mean(mydat < breaks[1])
  return(data.frame(x=breaks, y=c(y0, cumsum(h$counts)/length(mydat))))
  #return(data.frame(x=breaks, y=sapply(breaks, FUN=function(x) mean(mydat <= x))))
}
