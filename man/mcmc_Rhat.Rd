%% File Name: mcmc_Rhat.Rd
%% File Version: 0.07

\name{mcmc_Rhat}
\alias{mcmc_Rhat}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computation of the Rhat Statistic from a Single MCMC Chain
}
\description{
Computes the Rhat statistic from a single MCMC chain.
}
\usage{
mcmc_Rhat(mcmc_object, n_splits=3)
}
\arguments{
  \item{mcmc_object}{
Object of class \code{mcmc}
}
  \item{n_splits}{
Number of splits for MCMC chain
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Numeric vector
}

%\references{
%}

%\author{
%Alexander Robitzsch
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%}


\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Computation Rhat statistic for 2PNO model fitting by MCMC
#############################################################################

data(data.read)

# estimate 2PNO with MCMC with 3000 iterations and 500 burn-in iterations
mod <- sirt::mcmc.2pno( dat=data.read, iter=1000, burnin=100 )
# plot MCMC chains
plot( mod$mcmcobj, ask=TRUE )
# compute Rhat statistics
round( sirt::mcmc_Rhat( mod$mcmcobj[[1]] ), 3 )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{Testlets}

