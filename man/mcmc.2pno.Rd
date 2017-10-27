%% File Name: mcmc.2pno.Rd
%% File Version: 1.24

\name{mcmc.2pno}
\alias{mcmc.2pno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MCMC Estimation of the Two-Parameter Normal Ogive Item Response Model
}
\description{
This function estimates the Two-Parameter normal ogive item response model 
by MCMC sampling (Johnson & Albert, 1999, p. 195ff.).
}
\usage{
mcmc.2pno(dat, weights=NULL , burnin = 500, iter = 1000, N.sampvalues = 1000, 
      progress.iter = 50, save.theta = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame with dichotomous item responses
}
\item{weights}{
	An optional vector with student sample weights
		}
  \item{burnin}{
Number of burnin iterations
}
  \item{iter}{
Total number of iterations
}
  \item{N.sampvalues}{
Maximum number of sampled values to save
}
  \item{progress.iter}{
Display progress every \code{progress.iter}-th iteration. If no progress
display is wanted, then choose \code{progress.iter} larger than \code{iter}.
}
  \item{save.theta}{
Should theta values be saved?
}
}
\details{
The two-parameter normal ogive item response model with a probit link 
function is defined by
\deqn{ P(X_{pi} = 1 | \theta_p ) = \Phi ( a_i \theta_p - b_i )
		\quad ,  \quad \theta_p \sim N(0,1) }
Note that in this implementation non-informative priors for the item
parameters are chosen (Johnson & Albert, 1999, p. 195ff.).		
}
\value{
A list of class \code{mcmc.sirt} with following entries:
  \item{mcmcobj}{Object of class \code{mcmc.list}}
  \item{summary.mcmcobj}{Summary of the \code{mcmcobj} object. In this
  summary the Rhat statistic and the mode estimate MAP is included.
  The variable \code{PercSEratio} indicates the proportion of the Monte Carlo
  standard error in relation to the total standard deviation of the
  posterior distribution.}
  \item{burnin}{Number of burnin iterations}
  \item{iter}{Total number of iterations}
  \item{a.chain}{Sampled values of \eqn{a_i} parameters}
  \item{b.chain}{Sampled values of \eqn{b_i} parameters}
  \item{theta.chain}{Sampled values of \eqn{\theta_p} parameters}
  \item{deviance.chain}{Sampled values of Deviance values}  
  \item{EAP.rel}{EAP reliability}
  \item{person}{Data frame with EAP person parameter estimates for 
  	  \eqn{\theta_p} and their corresponding posterior standard
  	deviations}
  \item{dat}{Used data frame}
  \item{weights}{Used student weights}
  \item{\dots}{Further values}
}
% "EAP.rel" = EAP.rel , "dat" = dat0 , "weights" = weights , 

\references{
Johnson, V. E., & Albert, J. H. (1999). \emph{Ordinal Data Modeling}.
New York: Springer.
}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
S3 methods: \code{\link{summary.mcmc.sirt}}, \code{\link{plot.mcmc.sirt}}

For estimating the 2PL model with marginal maximum likelihood see
\code{\link{rasch.mml2}} or \code{\link{smirt}}.

A hierarchical version of this model can be estimated with
\code{\link{mcmc.2pnoh}}.
}
\examples{
\dontrun{	
#############################################################################
# EXAMPLE 1: Dataset Reading
#############################################################################
data(data.read)
# estimate 2PNO with MCMC with 3000 iterations and 500 burn-in iterations
mod <- sirt::mcmc.2pno( dat=data.read , iter=3000 , burnin=500 )
# plot MCMC chains
plot( mod$mcmcobj , ask=TRUE )
# write sampled chains into codafile
mcmclist2coda( mod$mcmcobj , name = "dataread_2pno" )
# summary
summary(mod)

#############################################################################
# EXAMPLE 2
#############################################################################
# simulate data
N <- 1000
I <- 10
b <- seq( -1.5 , 1.5 , len=I )
a <- rep( c(1,2) , I/2 )
theta1 <- stats::rnorm(N)
dat <- sirt::sim.raschtype( theta=theta1 , fixed.a =a , b=b )

#***
# Model 1: estimate model without weights
mod1 <- sirt::mcmc.2pno( dat , iter= 1500 , burnin=500)
mod1$summary.mcmcobj
plot( mod1$mcmcobj , ask=TRUE )

#***
# Model 2: estimate model with weights
# define weights
weights <- c( rep( 5 , N/4 ) , rep( .2 , 3/4*N ) )
mod2 <- sirt::mcmc.2pno( dat , weights=weights , iter= 1500 , burnin=500)
mod1$summary.mcmcobj
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Markov Chain Monte Carlo (MCMC)}
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
