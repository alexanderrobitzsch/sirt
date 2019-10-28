%% File Name: tetrachoric2.Rd
%% File Version: 1.363

\name{tetrachoric2}
\alias{tetrachoric2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tetrachoric Correlation Matrix
}
\description{
This function estimates a tetrachoric correlation matrix according to
the maximum likelihood estimation of Olsson (Olsson, 1979; \code{method="Ol"}),
the Tucker method (Method 2 of Froemel, 1971; \code{method="Tu"}) and
Divgi (1979, \code{method="Di"}).
In addition, an alternative non-iterative approximation of
Bonett and Price (2005; \code{method="Bo"}) is provided.
}
\usage{
tetrachoric2(dat, method="Ol", delta=0.007, maxit=1000000, cor.smooth=TRUE,
   progress=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
A data frame of dichotomous response
}
\item{method}{Computation method for calculating the tetrachoric correlation.
    The ML method is \code{method="Ol"} (which is the default),
    the Tucker method is \code{method="Tu"},
    the Divgi method is \code{method="Di"}
    the method of Bonett and
    Price (2005) is \code{method="Bo"}.}
  \item{delta}{
The step parameter. It is set by default to \eqn{2^{-7}}
which is approximately .007.
}
  \item{maxit}{
Maximum number of iterations.
}
\item{cor.smooth}{
Should smoothing of the tetrachoric correlation matrix
be performed to ensure positive definiteness? Choosing
\code{cor.smooth=TRUE}, the function \code{cor.smooth}
from the \pkg{psych} package is used for obtaining a positive definite
tetrachoric correlation matrix.
}
\item{progress}{Display progress? Default is \code{TRUE}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with following entries
\item{tau}{Item thresholds}
\item{rho}{Tetrachoric correlation matrix}
}


\references{
Bonett, D. G., & Price, R. M. (2005). Inferential methods for the tetrachoric
correlation coefficient. \emph{Journal of Educational and Behavioral Statistics,
30}(2), 213-225. \doi{10.3102/10769986030002213}

Divgi, D. R. (1979). Calculation of the tetrachoric correlation coefficient.
\emph{Psychometrika, 44}(2), 169-172.
\doi{10.1007/BF02293968}

Froemel, E. C. (1971). A comparison of computer routines for the
calculation of the tetrachoric correlation coefficient.
\emph{Psychometrika, 36}(2), 165-174.
\doi{10.1007/BF02291396}

Olsson, U. (1979). Maximum likelihood estimation of the polychoric correlation
coefficient. \emph{Psychometrika, 44}(4), 443-460.
\doi{10.1007/BF02296207}
}


\author{
Alexander Robitzsch

The code is adapted from an \R script of Cengiz Zopluoglu.
See \emph{http://sites.education.miami.edu/zopluoglu/software-programs/}.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also the \code{\link[psych:tetrachoric]{psych::tetrachoric}}
function in the \pkg{psych} package
 and the function \code{irtoys::tet} in the \pkg{irtoys} package.

See \code{\link{polychoric2}} for estimating polychoric correlations.
}
\examples{
#############################################################################
# EXAMPLE 1: data.read
#############################################################################

data(data.read)

# tetrachoric correlation from psych package
library(psych)
t0 <- psych::tetrachoric( data.read )$rho
# Olsson method (maximum likelihood estimation)
t1 <- sirt::tetrachoric2( data.read )$rho
# Divgi method
t2 <- sirt::tetrachoric2( data.read, method="Di"  )$rho
# Tucker method
t3 <- sirt::tetrachoric2( data.read, method="Tu" )$rho
# Bonett method
t4 <- sirt::tetrachoric2( data.read, method="Bo" )$rho

# maximum absolute deviation ML method
max( abs( t0 - t1 ) )
  ##   [1] 0.008224986
# mean absolute deviation Divgi method
max( abs( t0 - t2 ) )
  ##   [1] 0.1766688
# mean absolute deviation Tucker method
max( abs( t0 - t3 ) )
  ##   [1] 0.1766292
# mean absolute deviation Bonett method
max( abs( t0 - t4 ) )
  ##   [1] 0.05695522
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Tetrachoric correlation}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
