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
pairwise_association(mat, method = "infotheo::condentropy")
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (mat, method = "condentropy") 
{
    nsp <- nrow(mat)
    m_spxsp <- outer(1:nsp, 1:nsp, FUN = Vectorize(function(i, 
        j) {
        do.call(method, args = list(X = mat[i, ], Y = mat[j, 
            ]))
    }))
    diag(m_spxsp) <- NA
    return(m_spxsp)
  }
}