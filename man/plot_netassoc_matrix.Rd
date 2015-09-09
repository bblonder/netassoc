\name{plot_netassoc_matrix}
\alias{plot_netassoc_matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots matrix with colormap
}
\description{
Plots species x species or species x site matrix with color map
}
\usage{
plot_netassoc_matrix(data, 
  colors, onesided=FALSE, main="", 
  legend=TRUE, axis=TRUE, title=TRUE, cex.axis=0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Input matrix; assumed to have dimension names
}
  \item{colors}{
Vector of colors
}
  \item{onesided}{
If \code{TRUE}, assumes that \code{colors} maps to an entirely positive or entirely negative scale, i.e. without zero-crossings.
}
  \item{main}{
Title of plot.
}
  \item{legend}{
If \code{TRUE}, plots a colorbar legend.
}
  \item{axis}{
If \code{TRUE}, labels axes with dimension names.
}
  \item{title}{
If \code{TRUE}, plots a title.
}

  \item{cex.axis}{
Expansion factor for axis labels.
}
}
\value{
None; used for the side effect of making a plot.
}
\examples{
nsp <- 10
nsites <- 30
obs <- matrix(rpois(n=nsp*nsites,10),
  nrow=nsp,ncol=nsites,
  dimnames=list(paste("Species",1:nsp),paste("Site",1:nsites)))
plot_netassoc_matrix(obs, onesided=TRUE, col=heat.colors(5))

int <- matrix(rnorm(n=nsp^2),
  nrow=nsp,ncol=nsp,
  dimnames=list(paste("Species",1:nsp),paste("Species",1:nsp)))
plot_netassoc_matrix(int, onesided=FALSE, 
  col=colorRampPalette(c("red","white","blue"))(50))
}
