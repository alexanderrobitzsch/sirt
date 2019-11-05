%% File Name: polychoric2.Rd
%% File Version: 0.296

\name{polychoric2}
\alias{polychoric2}
\alias{sirt_rcpp_polychoric2}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Polychoric Correlation
}
\description{
This function estimates the polychoric correlation coefficient
using maximum likelihood estimation (Olsson, 1979).
}
\usage{
polychoric2(dat, maxiter=100, cor.smooth=TRUE, use_pbv=1, conv=1e-10,
      rho_init=NULL, weights=NULL)

## exported Rcpp function
sirt_rcpp_polychoric2( dat, maxK, maxiter, use_pbv, conv, rho_init, weights)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
A dataset with integer values \eqn{0,1,\ldots,K}
}
  \item{maxiter}{
Maximum number of iterations
}
  \item{cor.smooth}{
An optional logical indicating whether the polychoric correlation
matrix should be smooth to ensure positive definiteness.
}
\item{use_pbv}{Integer indicating whether the \pkg{pbv} package is used
for computation of bivariate normal distribution. \code{0} stands for
the simplest approximation in \pkg{sirt} (Cox & Wermuth, 1991, as
implemented in \code{\link{polychoric2}})
while versions \code{1} and
\code{2} uses the algorithm of \pkg{pbv} (the first one copied into
the \pkg{sirt} package, the second one linking Rcpp code to \pkg{pbv}.)}
\item{conv}{Convergence criterion}
\item{rho_init}{Optional matrix of initial values for polychoric correlations}
\item{weights}{Optional vector of sampling weights}
\item{maxK}{Maximum number of categories}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
A list with following entries
\item{tau}{Matrix of thresholds}
\item{rho}{Polychoric correlation matrix}
\item{Nobs}{Sample size for every item pair}
\item{maxcat}{Maximum number of categories per item}
}


\references{
Cox, D. R., & Wermuth, N. (1991). A simple approximation for bivariate and
trivariate normal integrals.
\emph{International Statistical Review, 59}(2), 263-269.

Olsson, U. (1979). Maximum likelihood estimation of the polychoric correlation
coefficient. \emph{Psychometrika, 44}(4), 443-460.
\doi{10.1007/BF02296207}
}
% Cox: \doi{10.2307/1403446}

%\author{
%Alexander Robitzsch
%}



%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See the \code{\link[psych:polychoric]{psych::polychoric}}
function in the \pkg{psych} package.

For estimating tetrachoric correlations see \code{\link{tetrachoric2}}.
}


\examples{
#############################################################################
# EXAMPLE 1: data.Students | activity scale
#############################################################################

data(data.Students, package="CDM")
dat <- data.Students[, paste0("act", 1:5 ) ]

# tetrachoric correlation from psych package
library(psych)
t0 <- psych::polychoric(dat)$rho
# Olsson method (maximum likelihood estimation)
t1 <- sirt::polychoric2(dat)$rho
# maximum absolute difference
max( abs( t0 - t1 ) )
  ##   [1] 0.004102429
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Polychoric correlation}
%% \keyword{Tetrachoric correlation}
