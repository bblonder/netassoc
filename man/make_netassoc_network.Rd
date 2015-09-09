\name{make_netassoc_network}
\alias{make_netassoc_network}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Infer species-association network
}
\description{
Infers a species association network by determining which co-occurrence patterns between species are more or less likely than expected under a null model of community assembly.   
}
\usage{
make_netassoc_network(obs, nul = vegan::permatfull(obs)$perm[[1]], 
  method="shrinkage", p.method="fdr", alpha=0.05, numnulls=1000, 
  plot=TRUE,plot.legend=TRUE, plot.title=TRUE, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{
A m x n community matrix describing the abundance of m species at n sites. Represents the observed data.
}
  \item{nul}{
A m x n community matrix describing the abundance of m species at n sites. Represents the regional null expectation data. The default value is a resampling of the observed data that preserves row and column sums.
}
  \item{method}{
The method used to calculate relationships between species. See \code{\link{partial_correlation}} for options. Defaults to a robust shrinkage estimator for inverse covariance matrices.
}
  \item{p.method}{
The method used to correct p-values for multiple comparisons. See \code{\link{p.adjust}} for options.
}
  \item{alpha}{
Analysis-wide Type I error rate, controlled via the argument \code{p.method}.
}
  \item{numnulls}{
Number of resamples of the \code{nul} matrix used to assemble null communities. Larger values produce more accurate results.
}
  \item{plot}{
If \code{TRUE}, plots all intermediate matrices calculated by the algorithm. Can be used to visualize input and output.
}
  \item{plot.title}{
If \code{TRUE}, adds titles to diagnostic plots.
}
  \item{plot.legend}{
If \code{TRUE}, adds legends to diagnostic plots.
}
  \item{verbose}{
If \code{TRUE}, prints status updates and progress bars during calculations.
}
}
\details{
Steps taken are:

1) obtaining input data and trimming to eliminate species that do not occur in any site
2) resampling a set of null community matrices from the expectation with the same richness and abundance as the observed community
3) calculating species co-occurrence scores for each pair of species within the observed matrix and all resampled null matrices
4) calculating standardized effect sizes and p-values for species' co-occurrence scores
5) thresholding effect sizes to retain only significant associations
6) converting matrix of scores to association network

The resulting network can be analyzed using functions from the \code{igraph} network package.

This process builds a Gaussian graphical model via estimating an inverse covariance matrix (precision matrix, which can be used to calculate partial correlation coefficients) for all species pairs. This graph is then compared to a distribution of null graphs, such that the final output is a graph with edge weights corresponding to standardized effect sizes after correction for multiple comparisons. 

A range of different methods are provided in \code{\link{partial_correlation}} for estimating relationships between species. Note that while a method is provided for the graphical lasso (L1-regularization) its use is not recommended, as it will produce very sparse null networks and then a narrow (or singular) distribution of null edge weights.
}
\value{
A list with the following components:
\item{\code{matrix_spsite_obs}}{Trimmed \code{obs} matrix}
\item{\code{matrix_spsite_nul}}{Trimmed \code{nul} matrix}
\item{\code{matrix_spsp_obs}}{Observed co-occurrence scores for all species}
\item{\code{matrix_spsp_ses_thresholded}}{Observed co-occurrence scores for all species after removing those with non-significant p-values}
\item{\code{matrix_spsp_pvalue}}{P-values for all species after correction for multiple comparisons}
\item{\code{network_all}}{An \code{igraph} object representing the association network}
\item{\code{network_pos}}{An \code{igraph} object representing an association network including only positive associations}
\item{\code{network_pos}}{An \code{igraph} object representing an association network including only negative associations}
}
\references{
Morueta-Holme, N., Blonder, B., et al. in preparation.
}
\seealso{vegan::permat}

\examples{
set.seed(1)
nsp <- 10
nsi <- 50
m_obs <- floor(matrix(rpois(nsp*nsi,lambda=5),ncol=nsi,nrow=nsp))
m_nul <- floor(matrix(rpois(nsp*nsi,lambda=5),ncol=nsi,nrow=nsp))

m_obs[1,1:(nsi/2)] <- rpois(n=nsi/2,lambda=20)
m_obs[2,1:(nsi/2)] <- rpois(n=nsi/2,lambda=20)

n <- make_netassoc_network(m_obs, m_nul,
  method="shrinkage",p.method='fdr', 
  numnulls=100, plot=TRUE,alpha=0.05)

n$network_all
}