# generate null sp x site matrices with same total abundance as observed and species probabilities drawn abundances within the base null matrix
generate_nul_resample <- function(nul, obs)
{
  nul_resample <- nul
  for (i in 1:ncol(nul_resample))
  {
    # zero out all values
    nul_resample[,i] <- 0
    
    # fill in resamples with the same total abundance
    samples <- sample(as.character(dimnames(nul)[[1]]),size=sum(obs[,i]),replace=TRUE,prob=nul[,i])
    samples <- factor(samples)
    abund_simulated <- tabulate(as.numeric(samples))
    nul_resample[as.character(levels(abund_simulated)),i] <- abund_simulated  	
  }	
  
  return(nul_resample)
}  

makenetwork <- function(obs, nul, whichmethod='pearson', kappa=2, numnulls=1000, plot=FALSE,verbose=TRUE)
{	
  if (is.matrix(obs) & is.null(dimnames(obs)))
  {
    dimnames(obs) <- list(1:nrow(obs),1:ncol(obs))
  }
  if (is.matrix(nul) & is.null(dimnames(nul)))
  {
    dimnames(nul) <- list(1:nrow(nul),1:ncol(nul))
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
  	stop(sprintf("Some sites are empty: %s",paste(names(problemsites),collapse=", ")))
  }

  if (plot==TRUE)
  {
    par(mfrow=c(2,4)) 
    par(cex.lab=0.5)
    par(cex.main=0.6)
    par(mar=c(0.5,2,2,0.5))
  }
  
  if (plot==TRUE)
  {
    plot_netassoc_matrix(obs, colors=colorRampPalette(c('white','green'))(51),onesided=TRUE,main="Observed sp x site")

    plot_netassoc_matrix(nul, colors=colorRampPalette(c('white','black'))(51),onesided=TRUE,main="Null sp x site")
  }
  
  
  if (verbose==TRUE) { cat('Generating null replicates...') }
  nulls <- vector(mode="list",length=numnulls)
  for (i in 1:numnulls)
  {
    if (verbose==TRUE) { cat(sprintf("nulrep%d ",i)) }
    nulls[[i]] <- generate_nul_resample(nul, obs)
  }
  if (verbose==TRUE) { cat('...done.\n') }
  
  if (plot==TRUE)
  {
    plot_netassoc_matrix(nulls[[sample(numnulls, 1)]], colors=colorRampPalette(c('white','black'))(51),onesided=TRUE,main="Example null resample sp x site")
  }	
  
  fm_obs <- matrix(NA,nrow=nrow(obs),ncol=nrow(obs))
  fm_nul_mean <- matrix(NA,nrow=nrow(obs),ncol=nrow(obs))
  fm_nul_sd <- matrix(NA,nrow=nrow(obs),ncol=nrow(obs))
  
  finalmatrix <- matrix(NA,nrow=nrow(obs),ncol=nrow(obs))
  
  if (verbose==TRUE) { cat('Calculating co-occurrence scores...') }
  count <- 0
  for (i in 1:nrow(obs))
  {
    for (j in i:nrow(obs))
    {
      if (i!=j)
      {
        count <- count + 1
        if (verbose==TRUE) { cat (sprintf('%d %d %.3f ', i, j, count/(nrow(obs)*(nrow(obs)-1)/2))) }
        
        if (verbose==TRUE) { cat('obs') }
        veci_obs <- as.numeric(obs[i,])
        vecj_obs <- as.numeric(obs[j,])
        cor_obs <- cor(x=veci_obs, y=vecj_obs, method=whichmethod)
        fm_obs[i,j] <- cor_obs
        if (verbose==TRUE) { cat('. ') }

        if (!is.na(cor_obs))
        {
          cor_nul <- rep(NA, numnulls)
          
          for (k in 1:numnulls)
          {
            if (verbose==TRUE) { cat(sprintf('nul%d',k)) }
            veci_nul <- as.numeric(nul[i,])
            vecj_nul <- as.numeric(nul[j,])
            cor_nul[k] <- cor(x=veci_nul, y=vecj_nul, method=whichmethod)
            if (verbose==TRUE) { cat('. ') }
          }
          
          fm_nul_mean[i,j] <- mean(cor_nul,na.rm=T)
          fm_nul_sd[i,j] <- sd(cor_nul,na.rm=T)    
        }
        else
        {
          if (verbose==TRUE) { cat('nul NONE') }
        }
        if (verbose==TRUE) { cat('\n') }
      }
    }
  }
  if (verbose==TRUE) { cat('...done.\n') }
  
  # calculate standard effect size
  if (verbose==TRUE) { cat('Calculating standardized effect sizes...') }
  finalmatrix <- (fm_obs - fm_nul_mean) / fm_nul_sd
  if (verbose==TRUE) { cat('...done.\n') }
  
  if (plot==TRUE)
  {
    plot_netassoc_matrix(fm_obs, xlab='Species',ylab='Species',colors=colorRampPalette(c('red','white','blue'))(51),main="Observed co-occurrence score for sp x sp")
    plot_netassoc_matrix(fm_nul_mean, xlab='Species',ylab='Species',colors=colorRampPalette(c('red','white','blue'))(51),main="Null mean co-occurrence score for sp x sp")
    plot_netassoc_matrix(fm_nul_sd,   xlab='Species',ylab='Species',colors=colorRampPalette(c('white','gray'))(51),onesided=TRUE,main="Null s.d. co-occurrence score for sp x sp")
  }
  
  dimnames(finalmatrix) <- list(row.names(obs),row.names(obs))
  
  # trim out low-value nodes	
  if (verbose==TRUE) { cat('Applying kappa threshold...') }
  finalmatrix_trimmed <- finalmatrix
  finalmatrix_trimmed[abs(finalmatrix_trimmed) < kappa] <- 0
  finalmatrix_trimmed[is.na(finalmatrix_trimmed)] <- 0
  finalmatrix_trimmed[is.infinite(finalmatrix_trimmed)] <- 0
  if (verbose==TRUE) { cat('...done.\n') }

  if (plot==TRUE)
  {
    plot_netassoc_matrix(finalmatrix_trimmed, xlab='Species',ylab='Species',colorRampPalette(c('red','white','blue'))(51),main="S.E.S. co-occurrence score for sp x sp")
  }
  
  # convert to network representation
  if (verbose==TRUE) { cat('Building network...') }
  network_all <- graph.adjacency(finalmatrix_trimmed,mode='upper',weighted=T)
  if (verbose==TRUE) { cat('...done.\n') }
  
  if (plot==TRUE)
  {
    plot_netassoc_network(network_all)
    title('Association network')
  }
  
  if (plot==TRUE)
  {
    par(mfrow=c(1,1))
  }
  

  return(list(
    matrix_spsite_obs=obs,
    matrix_spsite_null=nul,
    matrix_spsite_nulls=nulls,
    matrix_spsp_obs=fm_obs,
    matrix_spsp_null_mean=fm_nul_mean,
    matrix_spsp_null_sd=fm_nul_sd,
    matrix_spsp_ses_all=finalmatrix,
    matrix_spsp_ses_thresholded=finalmatrix_trimmed,
    network_all=network_all
    ))

}




# plot network
plot_netassoc_network <- function(network, layout = layout.auto(network), 
                                  vertex.label = V(network)$name, 
                                  vertex.color = NA, 
                                  vertex.shape = "none",
                                  vertex.label.color = "black", 
                                  vertex.label.family = "sans",
                                  edge.width = NULL, 
                                  edge.color = NULL, 
                                  edge.arrow.size = 0.05, 
                                  vertex.label.cex = 0.5, 
                                  ...)
{    
  if(is.null(edge.width))
  {
    if(length(E(network)$weight)==0)
    {
      edge.width=1
    }
    else
    {
      edge.width=sqrt(abs(E(network)$weight))
    }
  }
  
  if(is.null(edge.color))
  {
    if(length(E(network)$weight)==0)
    {
      edge.color <- 'black'  
      zlmin <- -1
      zlmax <- 1
    }
    else
    {
      edge.color <- ifelse(E(network)$weight > 0, rgb(0,0,1,0.8),rgb(1,0,0,0.8))
      zlmax <- max(abs(E(network)$weight),na.rm=T)
      zlmin = -1*zlmax
    }
  }
  
  plot(network,
       layout=layout,
       vertex.label=vertex.label,
       edge.color=edge.color,
       edge.width=edge.width,
       vertex.color=vertex.color,
       vertex.label.color=vertex.label.color,
       vertex.shape=vertex.shape,
       edge.arrow.size=edge.arrow.size,
       vertex.label.cex=vertex.label.cex,
       vertex.label.family=vertex.label.family,
       ...)
  
  colors=colorRampPalette(c('red','white','blue'))(51)
  legend('topleft',adj=c(0,0),legend=format(c(zlmin,zlmin/2+zlmax/2,zlmax),digits=2),fill=c(colors[1],colors[ceiling(length(colors)/2)],colors[length(colors)]),bg='white',cex=0.5)
}

plot_netassoc_matrix <- function(data, colors, ylab='Species',xlab='Site',onesided=FALSE,main="")
{
  zlmax <- max(abs(as.numeric(data)),na.rm=T)
  if (is.infinite(zlmax))
  {
    zlmax <- 1
  }
  if (onesided==TRUE)
  {
    zlmin = 0
  }
  else
  {
    zlmin = -1*zlmax
  }
  
  image(t(data),col=colors,axes=F,zlim=c(zlmin, zlmax),main=main)
  
  mtext(side=2,text=ylab,cex=0.5)
  mtext(side=3,text=xlab,cex=0.5)
  legend('topleft',adj=c(0,0),legend=format(c(zlmin,zlmin/2+zlmax/2,zlmax),digits=2),fill=c(colors[1],colors[ceiling(length(colors)/2)],colors[length(colors)]),bg='white',cex=0.5)
  box()
}




