%% File Name: polychoric2.Rd
%% File Version: 0.14

\name{polychoric2}
\alias{polychoric2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Polychoric Correlation
}
\description{
This function estimates the polychoric correlation coefficient
using maximum likelihood estimation (Olsson, 1979).
}
\usage{
polychoric2(dat, maxiter=100, cor.smooth=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
A dataset with integer values
}
  \item{maxiter}{
Maximum number of iterations
}
  \item{cor.smooth}{
An optional logical indicating whether the polychoric correlation
matrix should be smooth to ensure positive definiteness.
}
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
Olsson, U. (1979). Maximum likelihood estimation of the polychoric correlation
coefficient. \emph{Psychometrika}, \bold{44}, 443-460.
}


\author{
Alexander Robitzsch
}

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
t0 <- psych::polychoric( dat )$rho
# Olsson method (maximum likelihood estimation)
t1 <- sirt::polychoric2( dat )$rho
# maximum absolute difference
max( abs( t0 - t1 ) )
  ##   [1] 0.006914892
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Polychoric correlation}
\keyword{Tetrachoric correlation}
