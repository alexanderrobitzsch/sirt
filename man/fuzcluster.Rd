%% File Name: fuzcluster.Rd
%% File Version: 0.22

\name{fuzcluster}
\alias{fuzcluster}
\alias{summary.fuzcluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clustering for Continuous Fuzzy Data
}
\description{
This function performs clustering for continuous fuzzy data
for which membership functions are assumed to be Gaussian
(Denoeux, 2013). The mixture is also assumed to be Gaussian and
(conditionally cluster membership) independent.
}
\usage{
fuzcluster(dat_m, dat_s, K=2, nstarts=7, seed=NULL, maxiter=100,
     parmconv=0.001, fac.oldxsi=0.75, progress=TRUE)

\method{summary}{fuzcluster}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat_m}{
Centers for individual item specific membership functions
}
  \item{dat_s}{
Standard deviations for individual item specific membership functions
}
  \item{K}{
Number of latent classes
}
  \item{nstarts}{
Number of random starts. The default is 7 random starts.
}
  \item{seed}{
Simulation seed. If one value is provided, then only one
start is performed.
}
  \item{maxiter}{
Maximum number of iterations
}
  \item{parmconv}{
Maximum absolute change in parameters
}
\item{fac.oldxsi}{Convergence acceleration factor which
should take values between 0 and 1. The default is 0.75.}
  \item{progress}{
An optional logical indicating whether iteration progress
should be displayed.
}
\item{object}{Object of class \code{fuzcluster}}
\item{\dots}{Further arguments to be passed}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
A list with following entries

\item{deviance}{Deviance}
\item{iter}{Number of iterations}
\item{pi_est}{Estimated class probabilities}
\item{mu_est}{Cluster means}
\item{sd_est}{Cluster standard deviations}
\item{posterior}{Individual posterior distributions of cluster membership}
\item{seed}{Simulation seed for cluster solution}
\item{ic}{Information criteria}
}

\references{
Denoeux, T. (2013). Maximum likelihood estimation from uncertain data
in the belief function framework.
\emph{IEEE Transactions on Knowledge and Data Engineering, 25}, 119-130.
}


%\author{
%Alexander Robitzsch
%}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{fuzdiscr}} for estimating discrete distributions for
fuzzy data.

See the \pkg{fclust} package for fuzzy clustering.
}


\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: 2 classes and 3 items
#############################################################################

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# simulate data (2 classes and 3 items)
set.seed(876)
library(mvtnorm)
Ntot <- 1000  # number of subjects
# define SDs for simulating uncertainty
sd_uncertain <- c( .2, 1, 2 )

dat_m <- NULL   # data frame containing mean of membership function
dat_s <- NULL   # data frame containing SD of membership function

# *** Class 1
pi_class <- .6
Nclass <- Ntot * pi_class
mu <- c(3,1,0)
Sigma <- diag(3)
# simulate data
dat_m1 <- mvtnorm::rmvnorm( Nclass, mean=mu, sigma=Sigma )
dat_s1 <- matrix( stats::runif( Nclass * 3 ), nrow=Nclass )
for ( ii in 1:3){ dat_s1[,ii] <- dat_s1[,ii] * sd_uncertain[ii] }
dat_m <- rbind( dat_m, dat_m1 )
dat_s <- rbind( dat_s, dat_s1 )

# *** Class 2
pi_class <- .4
Nclass <- Ntot * pi_class
mu <- c(0,-2,0.4)
Sigma <- diag(c(0.5, 2, 2 ) )
# simulate data
dat_m1 <- mvtnorm::rmvnorm( Nclass, mean=mu, sigma=Sigma )
dat_s1 <- matrix( stats::runif( Nclass * 3 ), nrow=Nclass )
for ( ii in 1:3){ dat_s1[,ii] <- dat_s1[,ii] * sd_uncertain[ii] }
dat_m <- rbind( dat_m, dat_m1 )
dat_s <- rbind( dat_s, dat_s1 )
colnames(dat_s) <- colnames(dat_m) <- paste0("I", 1:3 )

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# estimation

#*** Model 1: Clustering with 8 random starts
res1 <- sirt::fuzcluster(K=2,dat_m, dat_s, nstarts=8, maxiter=25)
summary(res1)
  ##  Number of iterations=22 (Seed=5090 )
  ##  ---------------------------------------------------
  ##  Class probabilities (2 Classes)
  ##  [1] 0.4083 0.5917
  ##
  ##  Means
  ##           I1      I2     I3
  ##  [1,] 0.0595 -1.9070 0.4011
  ##  [2,] 3.0682  1.0233 0.0359
  ##
  ##  Standard deviations
  ##         [,1]   [,2]   [,3]
  ##  [1,] 0.7238 1.3712 1.2647
  ##  [2,] 0.9740 0.8500 0.7523

#*** Model 2: Clustering with one start with seed 4550
res2 <- sirt::fuzcluster(K=2,dat_m, dat_s, nstarts=1, seed=5090 )
summary(res2)

#*** Model 3: Clustering for crisp data
#             (assuming no uncertainty, i.e. dat_s=0)
res3 <- sirt::fuzcluster(K=2,dat_m, dat_s=0*dat_s, nstarts=30, maxiter=25)
summary(res3)
  ##  Class probabilities (2 Classes)
  ##  [1] 0.3645 0.6355
  ##
  ##  Means
  ##           I1      I2      I3
  ##  [1,] 0.0463 -1.9221  0.4481
  ##  [2,] 3.0527  1.0241 -0.0008
  ##
  ##  Standard deviations
  ##         [,1]   [,2]   [,3]
  ##  [1,] 0.7261 1.4541 1.4586
  ##  [2,] 0.9933 0.9592 0.9535

#*** Model 4: kmeans cluster analysis
res4 <- stats::kmeans( dat_m, centers=2 )
  ##   K-means clustering with 2 clusters of sizes 607, 393
  ##   Cluster means:
  ##             I1        I2          I3
  ##   1 3.01550780  1.035848 -0.01662275
  ##   2 0.03448309 -2.008209  0.48295067
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Fuzzy data}
\keyword{Clustering}
\keyword{summary}% __ONLY ONE__ keyword per line
