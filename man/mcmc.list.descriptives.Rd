%% File Name: mcmc.list.descriptives.Rd
%% File Version: 0.27

\name{mcmc.list.descriptives}
\alias{mcmc.list.descriptives}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computation of Descriptive Statistics for a \code{mcmc.list} Object
}
\description{
Computation of descriptive statistics, Rhat convergence statistic
and MAP for a \code{mcmc.list} object. The Rhat statistic
is computed by splitting one Monte Carlo chain into three segments of equal
length. The MAP is the mode estimate of the posterior distribution which is
approximated by the mode of the kernel density estimate.
}
\usage{
mcmc.list.descriptives( mcmcobj, quantiles=c(.025,.05,.1,.5,.9,.95,.975) )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmcobj}{
Object of class \code{mcmc.list}
}
\item{quantiles}{Quantiles to be calculated for all parameters}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A data frame with descriptive statistics for all parameters in
the \code{mcmc.list} object.
}
%\references{
%% ~put references to the literature/web site here ~
%}

%\author{
%Alexander Robitzsch
%}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{mcmclist2coda}} for writing an object of class \code{mcmc.list}
into a coda file (see also the \pkg{coda} package).
}

\examples{
\dontrun{
miceadds::library_install("coda")
miceadds::library_install("R2WinBUGS")

#############################################################################
# EXAMPLE 1: Logistic regression
#############################################################################

#***************************************
# (1) simulate data
set.seed(8765)
N <- 500
x1 <- stats::rnorm(N)
x2 <- stats::rnorm(N)
y <- 1*( stats::plogis( -.6 + .7*x1 + 1.1 *x2 ) > stats::runif(N) )

#***************************************
# (2) estimate logistic regression with glm
mod <- stats::glm( y ~ x1 + x2, family="binomial" )
summary(mod)

#***************************************
# (3) estimate model with rcppbugs package
b <- rcppbugs::mcmc.normal( stats::rnorm(3),mu=0,tau=0.0001)
y.hat <- rcppbugs::deterministic(function(x1,x2,b) {
             stats::plogis( b[1] + b[2]*x1 + b[3]*x2 ) }, x1, x2, b)
y.lik <- rcppbugs::mcmc.bernoulli( y, p=y.hat, observed=TRUE)
m <- rcppbugs::create.model(b, y.hat, y.lik)

#*** estimate model in rcppbugs; 5000 iterations, 1000 burnin iterations
ans <- rcppbugs::run.model(m, iterations=5000, burn=1000, adapt=1000, thin=5)
print(rcppbugs::get.ar(ans))     # get acceptance rate
print(apply(ans[["b"]],2,mean))  # get means of posterior

#*** convert rcppbugs into mcmclist object
mcmcobj <- data.frame( ans$b  )
colnames(mcmcobj) <- paste0("b",1:3)
mcmcobj <- as.matrix(mcmcobj)
class(mcmcobj) <- "mcmc"
attr(mcmcobj, "mcpar") <- c( 1, nrow(mcmcobj), 1 )
mcmcobj <- coda::as.mcmc.list( mcmcobj )

# plot results
plot(mcmcobj)

# summary
summ1 <-  sirt::mcmc.list.descriptives( mcmcobj )
summ1
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
