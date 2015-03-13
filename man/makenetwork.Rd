\name{makenetwork}
\alias{makenetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Infer species-association network
}
\description{
Infers a species association network by determining which co-occurrence patterns between species are more or less likely than expected under a null model of community assembly. Steps taken are:

1) obtaining input data and trimming to eliminate species that do not occur in any site
2) resampling a set of null community matrices from the expectation with the same richness and abundance as the observed community
3) calculating species co-occurrence scores for each pair of species within the observed matrix and all resampled null matrices
4) calculating standardized effect sizes for species' co-occurrence scores
5) thresholding effect sizes to retain only significant associations
6) converting matrix of scores to association network

The resulting network can be analyzed using functions from the \code{igraph} network package.
}
\usage{
makenetwork(obs, nul, 
  whichmethod = "pearson", kappa = 2, numnulls = 1000, 
  plot = FALSE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{
A m x n community matrix describing the abundance of m species at n sites. Represents the observed data.
}
  \item{nul}{
A m x n community matrix describing the abundance of m species at n sites. Represents the regional null expectation data.
}
  \item{whichmethod}{
Method used to calculate co-occurrence score between abundance vectors of species across sites. Can be either \code{'pearson'} or \code{'spearman'}, as per arguments in \code{\link{cor}}
}
  \item{kappa}{
Threshold value determining the minimum absolute standardized effect size to be retained in the network. Set to zero to retain all associations.
}
  \item{numnulls}{
Number of resamples of the \code{nul} matrix used to assemble null communities. Larger values produce more accurate results.
}
  \item{plot}{
If \code{TRUE}, plots all intermediate matrices calculated by the algorithm. Can be used to visualize input.
}
  \item{verbose}{
If \code{TRUE}, prints status updates and progress bars during calculations.
}
}
\details{
Full explanation of the algorithm is provided in the referenced paper.
}
\value{
A list with the following components:
\item{\code{matrix_spsite_obs}}{The trimmed \code{obs} matrix}
\item{\code{matrix_spsite_null}}{The trimmed \code{nul} matrix}
\item{\code{matrix_spsp_obs}}{The observed co-occurrence scores between species}
\item{\code{matrix_spsp_null_mean}}{The mean null co-occurrence scores between species}
\item{\code{matrix_spsp_null_sd}}{The s.d. null co-occurrence scores between species}
\item{\code{matrix_spsp_ses_all}}{The raw standardized effect size co-occurrence scores between species}
\item{\code{matrix_spsp_ses_thresholded}}{The thresholded standardized effect size co-occurrence scores between species}
\item{\code{network_all}}{An \code{igraph} object representing the association network}
\item{\code{network_pos}}{An \code{igraph} object representing an association network including only positive associations}
\item{\code{network_pos}}{An \code{igraph} object representing an association network including only negative associations}
}
\note{
Warnings may be generated when calculating co-occurrence scores using the \code{cor} function in cases where a pair of species co-occur identically at all sites. These warnings can be ignored.
}
\references{
Morueta-Holme, N., Blonder, B., et al. in preparation.
}

\examples{
# generate random data
set.seed(1)
nsp <- 10
nsi <- 5
m_obs <- floor(matrix(rgamma(nsp*nsi,shape=5),ncol=nsi,nrow=nsp))
m_nul <- floor(matrix(rexp(nsp*nsi,rate=0.05),ncol=nsi,nrow=nsp))

n <- makenetwork(m_obs, m_nul, numnulls=50, plot=TRUE)

# extract network
n$network_all

# example statistics
degree(n$network_all)
}
