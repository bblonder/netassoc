makenetwork <- function(obs, nul=vegan::permatfull(obs)$perm[[1]], method="shrinkage", p.method="fdr", alpha=0.5, numnulls=1000, plot=TRUE,verbose=TRUE)
{	
  if (is.matrix(obs) & is.null(dimnames(obs)))
  {
    dimnames(obs) <- list(paste("Species", 1:nrow(obs)),paste("Site", 1:ncol(obs)))
  }
  if (is.matrix(nul) & is.null(dimnames(nul)))
  {
    dimnames(nul) <- list(paste("Species", 1:nrow(nul)),paste("Site", 1:ncol(nul)))
  }
  
  obs <- as.matrix(obs,rownames.force=TRUE)
  nul <- as.matrix(nul,rownames.force=TRUE)
  
  if(!all(dim(obs)==dim(nul)))
  {
    stop("obs and nul must be same dimensionalities.")
  }

  problemspecies <- which(rowSums(obs)==0 | rowSums(nul)==0 | is.na(rowSums(nul)) | is.na(rowSums(obs)))
  if (length(problemspecies) > 0)
  {
	  stop(sprintf("Some species do not occur in any sites: %s",paste(names(problemspecies),collapse=", ")))
  }
  
  problemsites <- which(colSums(obs)==0 | colSums(nul)==0 | is.na(colSums(nul)) | is.na(colSums(obs)))
  if (length(problemsites) > 0)
  {
  	stop(sprintf("Some sites have no individuals of any species: %s",paste(names(problemsites),collapse=", ")))
  }
  
  nsp <- nrow(obs)

  if (plot==TRUE)
  {
    par(mfrow=c(3,3)) 
    par(cex.lab=0.5)
    par(cex.main=0.6)
    par(mar=c(0.5,2,2,0.5))
  }
  
  if (plot==TRUE)
  {
    plot_netassoc_matrix(obs, colors=colorRampPalette(c('white','green'))(51),onesided=TRUE,main="Observed sp x site")
    plot_netassoc_matrix(nul, colors=colorRampPalette(c('white','black'))(51),onesided=TRUE,main="Null sp x site")
   }
  
  
  
  if (verbose==TRUE) { cat('Calculating observed co-occurrence scores...\n') }
  pcor_obs <- getpcor(obs, method)
  
  if (plot==TRUE)
  {
     plot_netassoc_matrix(pcor_obs, colors=colorRampPalette(c('red','white','blue'))(51),onesided=FALSE,main="Obs sp x sp")
  }

  # create null distributions
  pcor_nul_all <- array(dim=c(nsp, nsp, numnulls))
  
  if (method=="glasso")
  {
    warning("Graphical lasso may not produce non-singular null distribution - SES and p-values may be unreasonable.\nProceed with caution.")
  }
  for (k in 1:numnulls)
  {
    # Generating null
    if (verbose==TRUE) { cat(sprintf('Generating null replicate %d...\n', k)) }
    
    nul_resampled <- generate_nul_resample(nul, obs)
    
    pcor_nul <- getpcor(nul_resampled, method)
    
    pcor_nul_all[,,k] <- pcor_nul
  }
  
  # calculate how often obs values are less than nul
  if (verbose==TRUE) { cat('Calculating standardized effect sizes...\n') }
  
  absobs_minus_absnul <- -1*sweep(abs(pcor_nul_all),c(1,2), abs(pcor_obs), "-")
  
  pcor_nul_mean <- apply(pcor_nul_all, c(1,2), mean)
  
  if (plot==TRUE)
  {
    plot_netassoc_matrix(pcor_nul_mean, colors=colorRampPalette(c('red','white','blue'))(51),onesided=FALSE,main="Nul mean sp x sp")
  }
  
  
  
  pcor_nul_sd <- apply(pcor_nul_all, c(1,2), sd)
  
  pcor_ses <- (pcor_obs - pcor_nul_mean) / pcor_nul_sd
  
  pcor_pvalues <- apply(absobs_minus_absnul, c(1,2), function(x) { mean(x<0) })
  
  # control false discovery rate
  if (verbose==TRUE) { cat('Adjusting p-values for multiple comparisons...\n') }
  pcor_pvalues_adjusted <- matrix(p.adjust(pcor_pvalues, method=p.method, n=((nsp)*(nsp-1)/2)),nrow=nsp,ncol=nsp)

  pcor_ses_trimmed <- pcor_ses
  pcor_ses_trimmed[pcor_pvalues_adjusted > alpha] <- NA
  
  pcor_ses_trimmed_pos <- pcor_ses_trimmed
  pcor_ses_trimmed_pos[pcor_ses_trimmed_pos <= 0] <- NA
  
  pcor_ses_trimmed_neg <- pcor_ses_trimmed
  pcor_ses_trimmed_neg[pcor_ses_trimmed_neg >= 0] <- NA
  pcor_ses_trimmed_neg <- -1 * pcor_ses_trimmed_neg
  
  if (plot==TRUE)
  {
    plot_netassoc_matrix(pcor_pvalues_adjusted, xlab='Species',ylab='Species',gray(1-1/1:100),main="P-values (corrected)",onesided=TRUE)
    plot_netassoc_matrix(pcor_ses, xlab='Species',ylab='Species',colorRampPalette(c('red','white','blue'))(51),main="S.E.S. co-occurrence score for sp x sp (raw)")    
    plot_netassoc_matrix(pcor_ses_trimmed, xlab='Species',ylab='Species',colorRampPalette(c('red','white','blue'))(51),main="S.E.S. co-occurrence score for sp x sp (trimmed)")
  }

  
  # convert to network representation
  if (verbose==TRUE) { cat('Building network...\n') }
  # overwrite NAs to zeros
  pcor_ses_trimmed_net <- pcor_ses_trimmed
  pcor_ses_trimmed_net[is.na(pcor_ses_trimmed_net)] <- 0
  pcor_ses_trimmed_pos_net <- pcor_ses_trimmed_pos
  pcor_ses_trimmed_pos_net[is.na(pcor_ses_trimmed_pos_net)] <- 0
  pcor_ses_trimmed_neg_net <- pcor_ses_trimmed_neg
  pcor_ses_trimmed_neg_net[is.na(pcor_ses_trimmed_neg_net)] <- 0
  
  network_all <- graph.adjacency(pcor_ses_trimmed_net,mode='upper',weighted=T)
  network_pos <- graph.adjacency(pcor_ses_trimmed_pos_net,mode='upper',weighted=T)
  network_neg <- graph.adjacency(pcor_ses_trimmed_neg_net,mode='upper',weighted=T)

  if (plot==TRUE)
  {
    plotnetassoc(network_all)
    title('Association network')
  }
  
  if (plot==TRUE)
  {
    par(mfrow=c(1,1))
  }
  

  return(list(
    matrix_spsite_obs=obs,
    matrix_spsite_nul=nul,
    matrix_spsp_obs=pcor_obs,
    matrix_spsp_ses_all=pcor_ses,
    matrix_spsp_ses_thresholded=pcor_ses_trimmed,
    matrix_spsp_pvalue=pcor_pvalues_adjusted,
    network_all=network_all,
    network_pos=network_pos,
    network_neg=network_neg
    ))

}






