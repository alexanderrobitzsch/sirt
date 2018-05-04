%% File Name: sia.sirt.Rd
%% File Version: 0.16

\name{sia.sirt}
\alias{sia.sirt}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Statistical Implicative Analysis (SIA)
}


\description{
This function is a simplified implementation of statistical implicative
analysis (Gras & Kuntz, 2008) which aims at deriving implications
\eqn{X_i \rightarrow X_j}. This means that solving item \eqn{i} implies
solving item \eqn{j}.
}

\usage{
sia.sirt(dat, significance = 0.85)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame with dichotomous item responses
}
  \item{significance}{
Minimum implicative probability for inclusion of an arrow in the graph.
The probability can be interpreted as a kind of significance level, i.e.
higher probabilities indicate more probable implications.
}
}

\details{
The test statistic for selection an implicative relation follows
Gras and Kuntz (2008).
Transitive arrows (implications) are removed from the graph.
If some implications are symmetric, then only the more probable
implication will be retained.
}


\value{
A list with following entries

  \item{adj.matrix}{Adjacency matrix of the graph. Transitive and
  symmetric implications (arrows) have been removed.}
  \item{adj.pot}{Adjacency matrix including all powers, i.e. all direct and
  indirect paths from item \eqn{i} to item \eqn{j}.}
  \item{adj.matrix.trans}{Adjacency matrix including transitive arrows.}
  \item{desc}{List with descriptive statistics of the graph.}
  \item{desc.item}{Descriptive statistics for each item.}
  \item{impl.int}{Implication intensity (probability) as the basis
  for deciding the significance of an arrow}
  \item{impl.t}{Corresponding \eqn{t} values of \code{impl.int}}
   \item{impl.significance}{Corresponding \eqn{p} values (significancies) 
   of \code{impl.int}}
  \item{conf.loev}{Confidence according to Loevinger (see Gras & Kuntz, 2008).
  This values are just conditional probabilities \eqn{P( X_j=1|X_i=1)}.}
  \item{graph.matr}{Matrix containing all arrows. Can be used
  for example for the \pkg{Rgraphviz} package.}
  \item{graph.edges}{Vector containing all edges of the graph, e.g. for
  the \pkg{Rgraphviz} package.}
  \item{igraph.matr}{Matrix containing all arrows for the \pkg{igraph}
      package.}
  \item{igraph.obj}{An object of the graph for the \pkg{igraph} package.}
  
}

%  \item{comp1 }{xx}
%  \item{comp2 }{xx}
 

\references{
Gras, R., & Kuntz, P. (2008). An overview of the statistical implicative analysis 
(SIA) development. In R. Gras, E. Suzuki, F. Guillet, & F. Spagnolo (Eds.).
\emph{Statistical Implicative Analysis} (pp. 11-40). 
Springer, Berlin Heidelberg.
}


\author{
Alexander Robitzsch
}


\note{
For an implementation of statistical implicative analysis in the 
C.H.I.C. (Classification Hierarchique, Implicative et Cohesitive)
software.

See \url{http://www.ardm.eu/contenu/logiciel-d-analyse-de-donnees-chic}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also the \pkg{IsingFit} package for calculating a graph for
dichotomous item responses using the Ising model.
}

%% ~~objects to See Also as \code{\link{help}}, ~~~

\examples{
#############################################################################
# EXAMPLE 1: SIA for data.read
#############################################################################

data(data.read)
dat <- data.read

res <- sirt::sia.sirt(dat , significance=.85 )

#*** plot results with igraph package
library(igraph)
plot( res$igraph.obj ) # , vertex.shape="rectangle" , vertex.size=30 )

\dontrun{
#*** plot results with qgraph package
miceadds::library_install(qgraph)
qgraph::qgraph( res$adj.matrix )

#*** plot results with Rgraphviz package
# Rgraphviz can only be obtained from Bioconductor
# If it should be downloaded, select TRUE for the following lines
if (FALSE){
     source("http://bioconductor.org/biocLite.R")
     biocLite("Rgraphviz")
            }
# define graph
grmatrix <- res$graph.matr
res.graph <- new("graphNEL", nodes= res$graph.edges , edgemode="directed")
# add edges
RR <- nrow(grmatrix)
for (rr in 1:RR){
    res.graph <- Rgraphviz::addEdge(grmatrix[rr,1], grmatrix[rr,2], res.graph , 1)
                    }                                                          
# define cex sizes and shapes
V <- length(res$graph.edges)
size2 <- rep(16,V)
shape2 <- rep("rectangle" , V )
names(shape2) <- names(size2) <- res$graph.edges
# plot graph
Rgraphviz::plot( res.graph, nodeAttrs =list("fontsize" = size2 , "shape" = shape2) )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Statistical implicative analysis}
\keyword{Graphical modeling}
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
