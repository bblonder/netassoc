\name{plot_netassoc_network}
\alias{plot_netassoc_network}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots species association network
}
\description{
Draws a network of species associations. By default edge widths are proportional to association strength and edge color reflects assocation type (blue, positive; red, negative).
}
\usage{
plot_netassoc_network(network, layout = layout.fruchterman.reingold(network), 
  vertex.label = V(network)$name, 
  vertex.color = NA, 
  vertex.shape = "none",
  vertex.label.color = "black", 
  vertex.label.family = "sans",
  edge.width = NULL, 
  edge.color = NULL, 
  edge.arrow.size = 0.05, 
  vertex.label.cex = 0.5, 
  legend = TRUE,
   ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{network}{
An \code{igraph} object corresponding to the association network
}
  \item{layout}{
Graphical layout. See \code{igraph::layout}.
}
  \item{vertex.label}{
String labels for species.
}
  \item{edge.width}{
Edge widths for links between species.
}
  \item{edge.color}{
Edge colors for links between species.
}
  \item{vertex.color}{
Vertex colors for species.
}
  \item{vertex.label.color}{
Vertex label colors for species.
}
  \item{vertex.shape}{
Vertex shape for species.
}
  \item{edge.arrow.size}{
Edge arrow size for links between species.
}
  \item{vertex.label.cex}{
Vertex label expansion factor for species.
}
  \item{vertex.label.family}{
Vertex shape font family for species.
}
  \item{legend}{
If \code{TRUE}, plots a scale legend.
}
  \item{\dots}{
Other arguments to be passed to \code{\link{plot.igraph}}.
}
}

\examples{
# generate random data
set.seed(5)
nsp <- 10
nsi <- 5
m_obs <- floor(matrix(rgamma(nsp*nsi,shape=5),ncol=nsi,nrow=nsp))
m_nul <- floor(matrix(rexp(nsp*nsi,rate=0.05),ncol=nsi,nrow=nsp))

n <- make_netassoc_network(m_obs, m_nul, numnulls=100, plot=TRUE,alpha=0.5)

# plot
plot_netassoc_network(n$network_all)

# plot using circular layout
plot_netassoc_network(n$network_all, layout=layout.circle(n$network_all))
}
