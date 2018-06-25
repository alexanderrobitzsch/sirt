%% File Name: reliability.nonlinearSEM.Rd
%% File Version: 0.18

\name{reliability.nonlinearSEM}
\alias{reliability.nonlinearSEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of Reliability for Confirmatory Factor Analyses
Based on Dichotomous Data
}
\description{
This function estimates a model based reliability using confirmatory factor
analysis (Green & Yang, 2009).
}
\usage{
reliability.nonlinearSEM(facloadings, thresh, resid.cov=NULL, cor.factors=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{facloadings}{
Matrix of factor loadings
}
  \item{thresh}{
Vector of thresholds
}
\item{resid.cov}{Matrix of residual covariances}
  \item{cor.factors}{
Optional matrix of covariances (correlations) between factors. The
default is a diagonal matrix with variances of 1.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list. The reliability is the list element \code{omega.rel}
}

\references{
Green, S. B., & Yang, Y. (2009). Reliability of summed item
scores using structural equation modeling: An alternative to
coefficient alpha. \emph{Psychometrika}, \bold{74}, 155-167.
}
\author{
Alexander Robitzsch
}
\note{
This function needs the \pkg{mvtnorm} package.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
This function is used in \code{\link{greenyang.reliability}}.
}
\examples{
#############################################################################
# EXAMPLE 1: Reading data set
#############################################################################
data(data.read)
dat <- data.read
I <- ncol(dat)

# define item clusters
itemcluster <- rep( 1:3, each=4)
error.corr <- diag(1,ncol(dat))
for ( ii in 1:3){
    ind.ii <- which( itemcluster==ii )
    error.corr[ ind.ii, ind.ii ] <- ii
        }
# estimate the model with error correlations
mod1 <- sirt::rasch.pml3( dat, error.corr=error.corr)
summary(mod1)

# extract item parameters
thresh <- - matrix( mod1$item$a * mod1$item$b, I, 1 )
A <- matrix( mod1$item$a * mod1$item$sigma, I, 1 )
# extract estimated correlation matrix
corM <- mod1$eps.corrM
# compute standardized factor loadings
facA <- 1 / sqrt( A^2 + 1 )
resvar <- 1 - facA^2
covM <- outer( sqrt(resvar[,1]), sqrt(resvar[,1] ) ) * corM
facloadings <- A *facA

# estimate reliability
rel1 <- sirt::reliability.nonlinearSEM( facloadings=facloadings, thresh=thresh,
           resid.cov=covM)
rel1$omega.rel
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Reliability}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
