%% File Name: rmvn.Rd
%% File Version: 0.08

\name{rmvn}
\alias{rmvn}
\alias{ruvn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of a Multivariate Normal Distribution with Exact Moments
}
\description{
Simulates a dataset from a multivariate or univariate normal distribution that
exactly fulfils the specified mean vector and the covariance matrix.
}
\usage{
# multivariate normal distribution
rmvn(N, mu, Sigma, exact=TRUE)

# univariate normal distribution
ruvn(N, mean=0, sd=1, exact=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
Sample size
}
  \item{mu}{
Mean vector
}
  \item{Sigma}{
Covariance matrix
}
  \item{exact}{
Logical indicating whether \code{mu} and \code{Sigma} should be exactly reproduced.
}
\item{mean}{Numeric value for mean}
\item{sd}{Numeric value for standard deviation}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
A dataframe or a vector
}

%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{mvtnorm::rmvnorm}, \code{mgcv::rmvn}
}


\examples{
#############################################################################
# EXAMPLE 1: Simulate multivariate normal data
#############################################################################

# define covariance matrix and mean vector
rho <- .8
Sigma <- matrix(rho,3,3)
diag(Sigma) <- 1
mu <- c(0,.5,1)

#* simulate data
set.seed(87)
dat <- sirt::rmvn(N=200, mu=mu, Sigma=Sigma)
#* check means and covariances
stats::cov.wt(dat, method="ML")

\dontrun{
#############################################################################
# EXAMPLE 2: Simulate univariate normal data
#############################################################################

#* simulate data
x <- sirt::ruvn(N=20, mean=.5, sd=1.2, exact=TRUE)
# check results
stats::var(x)
sirt:::sirt_var(x)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
