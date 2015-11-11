\name{pairwise_association}
\alias{pairwise_association}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pairwise associations
}
\description{
Computes pairwise associations between every row (species) in a species x site matrix. Note that usage of this function is advantageous when non-symmetric association metrics are desired, but the pairwise computation will prevent accounting for indirect effects between species. As such this function should be considered preliminary, and its use experimental.
}
\usage{
pairwise_association(mat, method = "condentropy")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
A m x n (species x site) matrix
}
  \item{method}{
The name of a function to call to calculate an association score. Must take two vector arguments (X,Y) and return a single numeric value. Default argument uses conditional information entropy statistic, although other functions (e.g. Jaccard similarity) are possible.
}
}
\value{
A n x n (species x species) matrix with NA diagonal values. May be non-symmetric depending on the method used.
}
\examples{
nsp <- 10
nsi <- 50
m_obs <- floor(matrix(rpois(nsp*nsi,lambda=5),ncol=nsi,nrow=nsp))
m_obs[1,1:(nsi/2)] <- rpois(n=nsi/2,lambda=20)

spxsp <- pairwise_association(m_obs, method="condentropy")
image(spxsp)
}