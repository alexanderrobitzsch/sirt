%% File Name: rasch.va.Rd
%% File Version: 0.11

\name{rasch.va}
\alias{rasch.va}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of the Rasch Model with Variational Approximation
}
\description{
This function estimates the Rasch model by the estimation method
of variational approximation (Rijmen & Vomlel, 2008).
}
\usage{
rasch.va(dat, globconv = 0.001, maxiter = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame with dichotomous item responses
}
  \item{globconv}{
Convergence criterion for item parameters
}
  \item{maxiter}{
Maximal number of iterations
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
A list with following entries:
  \item{sig}{Standard deviation of the trait}
  \item{item}{Data frame with item parameters}
  \item{xsi.ij}{Data frame with variational parameters \eqn{\xi_{ij}}}
  \item{mu.i}{Vector with individual means \eqn{\mu_i}}
\item{sigma2.i}{Vector with individual variances \eqn{\sigma_i^2}}  
%% ...
%%    res <- list("sig" = sqrt(sig2) , 
%%            "item"=item , "xsi.ij"=xsi.ij , "mu.i"=mu.i , "sigma2.i" = sigma2.i )
}
\references{
Rijmen, F., & Vomlel, J. (2008). Assessing the performance of
variational methods for mixed logistic regression models.
\emph{Journal of Statistical Computation and Simulation},
\bold{78}, 765-779.
}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
#############################################################################
# EXAMPLE 1: Rasch model
#############################################################################
set.seed(8706)
N <- 5000
I <- 20
dat <- sirt::sim.raschtype( stats::rnorm(N,sd=1.3) , b= seq(-2,2,len=I) )

# estimation via variational approximation
mod1 <- sirt::rasch.va(dat)

# estimation via marginal maximum likelihood
mod2 <- sirt::rasch.mml2(dat)

# estmation via joint maximum likelihood
mod3 <- sirt::rasch.jml(dat)

# compare sigma
round( c( mod1$sig , mod2$sd.trait ) , 3 )
## [1] 1.222 1.314

# compare b
round( cbind( mod1$item$b , mod2$item$b , mod3$item$itemdiff) , 3 )
##         [,1]   [,2]   [,3]
##  [1,] -1.898 -1.967 -2.090
##  [2,] -1.776 -1.841 -1.954
##  [3,] -1.561 -1.618 -1.715
##  [4,] -1.326 -1.375 -1.455
##  [5,] -1.121 -1.163 -1.228
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Variational approximation}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
