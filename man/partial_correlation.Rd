\name{partial_correlation}
\alias{partial_correlation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Partial correlation coefficients
}
\description{
Estimates the inverse covariance matrix then uses this matrix to calculate partial correlation coefficents.

Assumes that matrix rows correspond to different variables of interest.

The one exception is if \code{method="correlation"}; see below for details. 
}
\usage{
partial_correlation(mat, method, verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\itemize{}
  \item{mat}{
Input matrix.
}
  \item{method}{One of the following:
\code{"glasso"}: L1-regularized estimation via graphical lasso. The analysis is run for a range of penalty parameters; then a single inverse covariance matrix is returned whose penalty parameter yields the with the best Extended Bayesian Information Criterion value. Uses huge::huge and huge::huge.select.
\code{"shrinkage"}: Shrinkage estimate. The shrinkage parameter is estimated via an analytic formula. Uses corpcor::invcov.shrink.
\code{"exact"}: Directly calculates inverse covariance matrix. Results may be sensitive to highly collinear inputs. Use not recommended. Uses base::solve.
\code{"correlation"}: Directly calculates correlation matrix. Does not return partial correlation coefficients as in all other cases. Use not recommended. Uses stats::cor.
}
  \item{verbose}{
Binary flag determining whether diagnostic output is shown.
}
}

\value{
Returns a m x m upper triangular matrix of partial correlation coefficients from an input m x n matrix. 
}
\examples{
# load highly collinear economic data time series
data(longley) 
longley_ss <- t(longley[,c(1:5,7)]) # put data in correct input format

colors <- colorRampPalette(c("red","white","blue"))(10)
pc_shrinkage <- partial_correlation(longley_ss,method="shrinkage")

image(pc_shrinkage,zlim=c(-1,1),col=colors)
}