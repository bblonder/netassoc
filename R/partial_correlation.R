partial_correlation <- function(mat, method)
{ 

  # calculate inverse covariance matrix
  if (method=="glasso")
  {
    invcov <- as.matrix(huge::huge.select(huge(t(mat),method="glasso"),criterion="ebic")$opt.icov)
  }
  else if (method=="shrinkage")
  {
    # use auto-selected lambda shrinkage parameter
    invcov <- corpcor::invcov.shrink(t(mat))
  }
  else if (method=="ridge")
  {
    # use auto-selected lambda penalty parameter based on approximate leave one out cross validation
    invcov <- rags2ridges::optPenalty.aLOOCV(t(mat), lambdaMin=1e-3,lambdaMax=1e4,step=100,graph=FALSE,verbose=FALSE)$optPrec
  }
  else if (method=="exact")
  {
    # estimate inverse covariance from 
    cat('Calculating exact inverse covariance...\n')
    invcov <- solve(cov(t(mat)))
  }
  
  # convert inverse covariances to partial correlation coefficients
  if (method!="correlation")
  {
    pcor <- decompose.invcov(invcov)$pr
  }
  else
  {
    pcor <- cor(t(mat))
  }
  pcor[lower.tri(pcor,diag=TRUE)] <- NA
  
  return(pcor)
}