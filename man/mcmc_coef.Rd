\name{mcmc_coef}
\alias{mcmc_coef}
\alias{mcmc_vcov}
\alias{mcmc_confint}
\alias{mcmc_plot}
\alias{mcmc_summary}
\alias{mcmc_derivedPars}
\alias{mcmc_WaldTest}
\alias{summary.mcmc_WaldTest}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Some Methods for Objects of Class \code{mcmc.list}
}

\description{
Some methods for objects of class \code{mcmc.list} created
from the \pkg{coda} package.
}

\usage{
## coefficients	
mcmc_coef(mcmcobj, exclude = "deviance")

## covariance matrix 
mcmc_vcov(mcmcobj, exclude = "deviance")

## confidence interval
mcmc_confint( mcmcobj, parm, level = .95, exclude="deviance" )

## summary function
mcmc_summary( mcmcobj , quantiles=c(.025,.05,.50,.95,.975) )

## plot function
mcmc_plot(mcmcobj, ...)

## inclusion of derived parameters in mcmc object
mcmc_derivedPars( mcmcobj , derivedPars )

## Wald test for parameters
mcmc_WaldTest( mcmcobj , hypotheses )

\method{summary}{mcmc_WaldTest}(object, digits=3, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmcobj}{
Objects of class \code{mcmc.list} as created by
\code{\link[coda:mcmc]{coda::mcmc}}
}
  \item{exclude}{
Vector of parameters which should be excluded in calculations
}
\item{parm}{Optional vector of parameters}
\item{level}{Confidence level}
\item{quantiles}{Vector of quantiles to be computed.}
\item{\dots}{Parameters to be passed to \code{mcmc_plot}.
See \code{LAM::plot.amh} for arguments.
}
\item{derivedPars}{List with derived parameters (see examples).}
\item{hypotheses}{List with hypotheses of the form
\eqn{g_i( \bold{\theta})=0}.}
\item{object}{Object of class \code{mcmc_WaldTest}.}
\item{digits}{Number of digits used for rounding.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

%\value{
%}
%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Alexander Robitzsch
}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[coda:mcmc]{coda::mcmc}}
}

\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Logistic regression in rcppbugs package 
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
mod <- stats::glm( y ~ x1 + x2 , family="binomial" )
summary(mod)

#***************************************
# (3) estimate model with rcppbugs package
library(rcppbugs)
b <- rcppbugs::mcmc.normal( stats::rnorm(3),mu=0,tau=0.0001)
y.hat <- rcppbugs::deterministic( function(x1,x2,b){
                stats::plogis( b[1] + b[2]*x1 + b[3]*x2 ) }, 
                  x1 , x2 , b)
y.lik <- rcppbugs::mcmc.bernoulli( y , p = y.hat, observed = TRUE)
model <- rcppbugs::create.model(b, y.hat, y.lik)

#*** estimate model in rcppbugs; 5000 iterations, 1000 burnin iterations
n.burnin <- 500 ; n.iter <- 2000 ; thin <- 2
ans <- rcppbugs::run.model(model , iterations=n.iter, burn=n.burnin, adapt=200, thin=thin)
print(rcppbugs::get.ar(ans)) # get acceptance rate
print(apply(ans[["b"]],2,mean)) # get means of posterior

#*** convert rcppbugs into mcmclist object
mcmcobj <- data.frame( ans$b )
colnames(mcmcobj) <- paste0("b",1:3)
mcmcobj <- as.matrix(mcmcobj)
class(mcmcobj) <- "mcmc"
attr(mcmcobj, "mcpar") <- c( n.burnin+1 , n.iter , thin )
mcmcobj <- coda::mcmc( mcmcobj )

# coefficients, variance covariance matrix and confidence interval
mcmc_coef(mcmcobj)
mcmc_vcov(mcmcobj)
mcmc_confint( mcmcobj , level = .90 )

# summary and plot
mcmc_summary(mcmcobj)
mcmc_plot(mcmcobj, ask=TRUE)

# include derived parameters in mcmc object
derivedPars <- list( "diff12" = ~ I(b2-b1) , "diff13" = ~ I(b3-b1) )
mcmcobj2 <- sirt::mcmc_derivedPars(mcmcobj , derivedPars = derivedPars )
mcmc_summary(mcmcobj2)

#*** Wald test for parameters
 # hyp1: b2 - 0.5 = 0
 # hyp2: b2 * b3 = 0
hypotheses <- list( "hyp1" = ~ I( b2 - .5 )  , "hyp2" = ~ I( b2*b3 ) )
test1 <- sirt::mcmc_WaldTest( mcmcobj , hypotheses=hypotheses )
summary(test1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
