%% File Name: mcmc.2pnoh.Rd
%% File Version: 0.19

\name{mcmc.2pnoh}
\alias{mcmc.2pnoh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MCMC Estimation of the Hierarchical IRT Model for Criterion-Referenced
Measurement
}
\description{
This function estimates the hierarchical IRT model for criterion-referenced
measurement which is based on a two-parameter normal ogive response
function (Janssen, Tuerlinckx, Meulders & de Boeck, 2000).
}
\usage{
mcmc.2pnoh(dat, itemgroups , prob.mastery=c(.5,.8) , weights=NULL , 
      burnin = 500, iter = 1000, N.sampvalues = 1000, 
      progress.iter = 50, prior.variance=c(1,1) , save.theta = FALSE)
}

\arguments{
  \item{dat}{
Data frame with dichotomous item responses
}
\item{itemgroups}{Vector with characters or integers which define the 
	criterion to which an item is associated.}
\item{prob.mastery}{Probability levels which define nonmastery, transition
	and mastery stage (see Details) }	
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
\item{prior.variance}{Scale parameter of the inverse gamma distribution
	for the \eqn{\sigma^2} and \eqn{\nu^2} item variance parameters}
  \item{save.theta}{
Should theta values be saved?
}
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETAILS
\details{
The hierarchical IRT model for criterion-referenced measurement
(Janssen et al., 2000) assumes that every item \eqn{i} intends
to measure a criterion \eqn{k}. The item response function is defined as
\deqn{ P(X_{pik} = 1 | \theta_p ) = 
		\Phi [ \alpha_{ik} ( \theta_p - \beta_{ik} ) ]
		\quad ,  \quad \theta_p \sim N(0,1) }
Item parameters \eqn{(\alpha_{ik},\beta_{ik})} are hierarchically modeled, i.e.
\deqn{ \beta_{ik} \sim N( \xi_k , \sigma^2 ) \quad \mbox{and} \quad 
	\alpha_{ik} \sim N( \omega_k , \nu^2 ) }
In the \code{mcmc.list} output object, also the derived parameters 
\eqn{d_{ik}= \alpha_{ik} \beta_{ik}} and \eqn{\tau_k = \xi_k \omega_k} are
calculated.
Mastery and nonmastery probabilities are based on a reference item \eqn{Y_{k}}
of criterion \eqn{k} and a response function
\deqn{ P(Y_{pk} = 1 | \theta_p ) = 
		\Phi [ \omega_{k} ( \theta_p - \xi_{k} ) ]
		\quad ,  \quad \theta_p \sim N(0,1) }
With known item parameters and person parameters, response probabilities of
criterion \eqn{k} are calculated. If a response probability of criterion \eqn{k}
is larger than \code{prob.mastery[2]}, then a student is defined as a 
master. If this probability is smaller than \code{prob.mastery[1]}, then
a student is a nonmaster. In all other cases, students are in a transition
stage.

In the \code{mcmcobj} output object, the parameters \code{d[i]} are defined by 
\eqn{d_{ik} = \alpha_{ik} \cdot \beta_{ik}} while \code{tau[k]} are defined by
\eqn{ \tau_k = \xi_k \cdot \omega_k }.
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VALUE

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
  \item{alpha.chain}{Sampled values of \eqn{\alpha_{ik}} parameters}
  \item{beta.chain}{Sampled values of \eqn{\beta_{ik}} parameters}
  \item{xi.chain}{Sampled values of \eqn{\xi_{k}} parameters}  
  \item{omega.chain}{Sampled values of \eqn{\omega_{k}} parameters}
  \item{sigma.chain}{Sampled values of \eqn{\sigma} parameter}
  \item{nu.chain}{Sampled values of \eqn{\nu} parameter}
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REFERENCES
\references{
Janssen, R., Tuerlinckx, F., Meulders, M., & de Boeck, P. (2000).
A hierarchical IRT model for criterion-referenced measurement.
\emph{Journal of Educational and Behavioral Statistics}, \bold{25}, 285-306.		
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

The two-parameter normal ogive model can be estimated with
\code{\link{mcmc.2pno}}.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMPLES

\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Simulated data according to Janssen et al. (2000, Table 2)
#############################################################################

N <- 1000
Ik <- c(4,6,8,5,9,6,8,6,5)
xi.k <- c( -.89 , -1.13 , -1.23 , .06 , -1.41 , -.66 , -1.09 , .57 , -2.44)
omega.k <- c(.98 , .91 , .76 , .74 , .71 , .80 , .79 , .82 , .54)

# select 4 attributes
K <- 4
Ik <- Ik[1:K] ; xi.k <- xi.k[1:K] ; omega.k <- omega.k[1:K]
sig2 <- 3.02
nu2 <- .09
I <- sum(Ik)
b <- rep( xi.k , Ik ) + stats::rnorm(I , sd = sqrt(sig2) )
a <- rep( omega.k , Ik ) + stats::rnorm(I , sd = sqrt(nu2) )
theta1 <- stats::rnorm(N)
t1 <- rep(1,N)
p1 <- stats::pnorm( outer(t1,a) * ( theta1 - outer(t1,b) ) )
dat <- 1  * ( p1 > stats::runif(N*I)  )
itemgroups <- rep( paste0("A" , 1:K ) , Ik )

# estimate model
mod <- sirt::mcmc.2pnoh(dat  , itemgroups , burnin=200 , iter=1000 )
# summary
summary(mod)
# plot
plot(mod$mcmcobj , ask=TRUE)
# write coda files
mcmclist2coda( mod$mcmcobj , name = "simul_2pnoh" )
}	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Markov Chain Monte Carlo (MCMC)}
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
