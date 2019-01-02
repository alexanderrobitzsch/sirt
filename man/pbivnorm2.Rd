%% File Name: pbivnorm2.Rd
%% File Version: 0.25

\name{pbivnorm2}
\alias{pbivnorm2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cumulative Function for the Bivariate Normal Distribution
}
\description{
This function evaluates the bivariate normal distribution
\eqn{\Phi_2 ( x, y ; \rho )}
assuming zero means and unit variances. It uses a simple approximation
by Cox and Wermuth (1991) with corrected formulas in Hong (1999).
}
\usage{
pbivnorm2(x, y, rho)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vector of \eqn{x} coordinates
}
  \item{y}{
Vector of \eqn{y} coordinates
}
  \item{rho}{
Vector of correlations between random normal variates
}
}
%% \details{
%%  ~~ If necessary, more details than the description above ~~
%% }

\value{
Vector of probabilities
}

\references{
Cox, D. R., & Wermuth, N. (1991). A simple approximation for bivariate and
trivariate normal integrals.
\emph{International Statistical Review, 59}, 263-269.

Hong, H. P. (1999). An approximation to bivariate and trivariate normal
integrals. \emph{Engineering and Environmental Systems, 16}, 115-127.
}



%\author{
%Alexander Robitzsch
%}


\note{
The function is less precise for correlations near 1 or -1.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also the
\code{\link[pbivnorm:pbivnorm]{pbivnorm::pbivnorm}}
function in the \pkg{pbivnorm} package.
}
\examples{
library(pbivnorm)
# define input
x <- c(0, 0,  .5, 1, 1  )
y <- c( 0, -.5,  1, 3, .5 )
rho <- c( .2, .8, -.4, .6, .5 )
# compare pbivnorm2 and pbivnorm functions
pbiv2 <- sirt::pbivnorm2( x=x, y=y, rho=rho )
pbiv <- pbivnorm::pbivnorm(  x,  y, rho=rho )
max( abs(pbiv-pbiv2))
  ## [1] 0.0030626
round( cbind( x, y, rho,pbiv, pbiv2 ), 4 )
  ##          x    y  rho   pbiv  pbiv2
  ##   [1,] 0.0  0.0  0.2 0.2820 0.2821
  ##   [2,] 0.0 -0.5  0.8 0.2778 0.2747
  ##   [3,] 0.5  1.0 -0.4 0.5514 0.5514
  ##   [4,] 1.0  3.0  0.6 0.8412 0.8412
  ##   [5,] 1.0  0.5  0.5 0.6303 0.6304
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Bivariate normal distribution}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
