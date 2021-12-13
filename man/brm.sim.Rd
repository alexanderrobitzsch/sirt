%% File Name: brm.sim.Rd
%% File Version: 0.533

\name{brm-Methods}
\alias{brm.sim}
\alias{brm.irf}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Functions for the Beta Item Response Model
}


\description{
Functions for simulating and estimating the Beta item response model
(Noel & Dauvier, 2007). \code{brm.sim} can be used for
simulating the model, \code{brm.irf} computes the item response
function. The Beta item response model is estimated as a discrete
version to enable estimation in \emph{standard} IRT software like
\pkg{mirt} or \pkg{TAM} packages.
}

\usage{
# simulating the beta item response model
brm.sim(theta, delta, tau, K=NULL)

# computing the item response function of the beta item response model
brm.irf( Theta, delta, tau, ncat, thdim=1, eps=1E-10 )
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
Ability vector of \eqn{\theta} values
}
  \item{delta}{
Vector of item difficulty parameters
}
  \item{tau}{
Vector item dispersion parameters
}
  \item{K}{
Number of discretized categories. The default is \code{NULL} which
means that the simulated item responses are real number values
between 0 and 1. If an integer \code{K} chosen, then values
are discretized such that values of 0, 1, ..., \eqn{K}-1 arise.
}
\item{Theta}{Matrix of the ability vector \eqn{\bold{\theta}}}
\item{ncat}{Number of categories}
\item{thdim}{Theta dimension in the matrix \code{Theta} on
which the item loads.}
\item{eps}{Nuisance parameter which stabilize probabilities.}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% DETAILS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\details{
The discrete version of the beta item response model is defined as follows.
Assume that for item \eqn{i} there are \eqn{K} categories resulting in
values \eqn{k=0,1,\dots,K-1}. Each value \eqn{k} is associated with a
corresponding the transformed value in \eqn{[0,1]}, namely
\eqn{ q (k)=1/(2 \cdot K), 1/(2 \cdot K) + 1/K, \ldots,  1 - 1/(2 \cdot K) }.
The item response model is defined as
\deqn{ P( X_{pi}=x_{pi} | \theta_p)  \propto
      q( x_{pi} )^{ m_{pi} - 1 } [ 1- q( x_{pi} ) ]^{ n_{pi} - 1 } }
This density is a discrete version of a Beta distribution with
shape parameters \eqn{m_{pi}} and \eqn{n_{pi}}. These parameters are
defined as
\deqn{ m_{pi}=\mathrm{exp} \left[ ( \theta_p - \delta_i + \tau_i ) / 2 \right]
    \qquad \mbox{and} \qquad
    n_{pi}=\mathrm{exp} \left[ ( - \theta_p + \delta_i + \tau_i ) / 2 \right]
                }

The item response function can also be formulated as
\deqn{ \mathrm{log} \left[ P( X_{pi}=x_{pi} | \theta_p) \right]  \propto
      ( m_{pi} - 1 ) \cdot \mathrm{log} [ q( x_{pi} ) ] +
      ( n_{pi} - 1 )  \cdot \mathrm{log} [ 1- q( x_{pi} ) ]
                      }

The item parameters can be reparameterized as
\eqn{ a_{i}=\mathrm{exp} \left[ ( - \delta_i + \tau_i ) / 2 \right]} and
\eqn{ b_{i}=\mathrm{exp} \left[ ( \delta_i + \tau_i ) / 2 \right]}.

Then, the original item parameters can be retrieved by
\eqn{\tau_i=\mathrm{log} ( a_i b_i)} and
\eqn{\delta_i=\mathrm{log} ( b_i / a_i)}.
Using \eqn{ \gamma _p=\mathrm{exp} ( \theta_p / 2) }, we obtain

\deqn{ \mathrm{log} \left[ P( X_{pi}=x_{pi} | \theta_p) \right]  \propto
       a_{i} \gamma_p  \cdot \mathrm{log} [ q( x_{pi} ) ] +
       b_i / \gamma_p   \cdot \mathrm{log} [ 1- q( x_{pi} ) ] -
      \left[ \mathrm{log} q( x_{pi} ) + \mathrm{log} [ 1- q( x_{pi} ) ] \right]
                      }

This formulation enables the specification of the Beta item response
model as a structured latent class model
(see \code{\link[TAM:tam.mml.3pl]{TAM::tam.mml.3pl}};
Example 1).

See Smithson and Verkuilen (2006) for motivations for treating
continuous indicators not as normally distributed variables.
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% VALUE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\value{
A simulated dataset of item responses if \code{brm.sim} is applied.

A matrix of item response probabilities if \code{brm.irf} is applied.
}


\references{
Gruen, B., Kosmidis, I., & Zeileis, A. (2012). Extended Beta regression
in \R: Shaken, stirred, mixed, and partitioned.
\emph{Journal of Statistical Software, 48}(11), 1-25.
\doi{10.18637/jss.v048.i11}

Noel, Y., & Dauvier, B. (2007). A beta item response model for continuous
bounded responses. \emph{Applied Psychological Measurement,
31}(1), 47-73. \doi{10.1177/0146621605287691}

Smithson, M., & Verkuilen, J. (2006). A better lemon squeezer?
Maximum-likelihood regression with beta-distributed dependent variables.
\emph{Psychological Methods, 11}(1), 54-71.
doi: 10.1037/1082-989X.11.1.54
}

%\author{
%Alexander Robitzsch
%}

%\note{
%
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\seealso{
See also the \pkg{betareg} package for fitting Beta regression
regression models in \R (Gruen, Kosmidis & Zeileis, 2012).
}

\examples{
#############################################################################
# EXAMPLE 1: Simulated data beta response model
#############################################################################

#*** (1) Simulation of the beta response model
# Table 3 (p. 65) of Noel and Dauvier (2007)
delta <- c( -.942, -.649, -.603, -.398, -.379, .523, .649, .781, .907 )
tau <- c( .382, .166, 1.799, .615, 2.092, 1.988, 1.899, 1.439, 1.057 )
K <- 5        # number of categories for discretization
N <- 500        # number of persons
I <- length(delta) # number of items

set.seed(865)
theta <- stats::rnorm( N )
dat <- sirt::brm.sim( theta=theta, delta=delta, tau=tau, K=K)
psych::describe(dat)

#*** (2) some preliminaries for estimation of the model in mirt
#*** define a mirt function
library(mirt)
Theta <- matrix( seq( -4, 4, len=21), ncol=1 )

# compute item response function
ii <- 1     # item ii=1
b1 <- sirt::brm.irf( Theta=Theta, delta=delta[ii], tau=tau[ii],  ncat=K )
# plot item response functions
graphics::matplot( Theta[,1], b1, type="l" )

#*** defining the beta item response function for estimation in mirt
par <- c( 0, 1,  1)
names(par) <- c( "delta", "tau","thdim")
est <- c( TRUE, TRUE, FALSE )
names(est) <- names(par)
brm.icc <- function( par, Theta, ncat ){
     delta <- par[1]
     tau <- par[2]
     thdim <- par[3]
     probs <- sirt::brm.irf( Theta=Theta, delta=delta, tau=tau,  ncat=ncat,
            thdim=thdim)
     return(probs)
            }
name <- "brm"
# create item response function
brm.itemfct <- mirt::createItem(name, par=par, est=est, P=brm.icc)
#*** define model in mirt
mirtmodel <- mirt::mirt.model("
           F1=1-9
            " )
itemtype <- rep("brm", I )
customItems <- list("brm"=brm.itemfct)

# define parameters to be estimated
mod1.pars <- mirt::mirt(dat, mirtmodel, itemtype=itemtype,
                   customItems=customItems, pars="values")

\dontrun{
#*** (3) estimate beta item response model in mirt
mod1 <- mirt::mirt(dat,mirtmodel, itemtype=itemtype, customItems=customItems,
               pars=mod1.pars, verbose=TRUE  )
# model summaries
print(mod1)
summary(mod1)
coef(mod1)
# estimated coefficients and comparison with simulated data
cbind( sirt::mirt.wrapper.coef( mod1 )$coef, delta, tau )
mirt.wrapper.itemplot(mod1,ask=TRUE)

#---------------------------
# estimate beta item response model in TAM
library(TAM)

# define the skill space: standard normal distribution
TP <- 21                   # number of theta points
theta.k <- diag(TP)
theta.vec <-  seq( -6,6, len=TP)
d1 <- stats::dnorm(theta.vec)
d1 <- d1 / sum(d1)
delta.designmatrix <- matrix( log(d1), ncol=1 )
delta.fixed <- cbind( 1, 1, 1 )

# define design matrix E
E <- array(0, dim=c(I,K,TP,2*I + 1) )
dimnames(E)[[1]] <- items <- colnames(dat)
dimnames(E)[[4]] <- c( paste0( rep( items, each=2 ),
        rep( c("_a","_b" ), I) ), "one" )
for (ii in 1:I){
    for (kk in 1:K){
      for (tt in 1:TP){
        qk <- (2*(kk-1)+1)/(2*K)
        gammap <- exp( theta.vec[tt] / 2 )
        E[ii, kk, tt, 2*(ii-1) + 1 ] <- gammap * log( qk )
        E[ii, kk, tt, 2*(ii-1) + 2 ] <- 1 / gammap * log( 1 - qk )
        E[ii, kk, tt, 2*I+1 ] <- - log(qk) - log( 1 - qk )
                    }
            }
        }
gammaslope.fixed <- cbind( 2*I+1, 1 )
gammaslope <- exp( rep(0,2*I+1) )

# estimate model in TAM
mod2 <- TAM::tam.mml.3pl(resp=dat, E=E,control=list(maxiter=100),
              skillspace="discrete", delta.designmatrix=delta.designmatrix,
              delta.fixed=delta.fixed, theta.k=theta.k, gammaslope=gammaslope,
              gammaslope.fixed=gammaslope.fixed, notA=TRUE )
summary(mod2)

# extract original tau and delta parameters
m1 <- matrix( mod2$gammaslope[1:(2*I) ], ncol=2, byrow=TRUE )
m1 <- as.data.frame(m1)
colnames(m1) <- c("a","b")
m1$delta.TAM <- log( m1$b / m1$a)
m1$tau.TAM <- log( m1$a * m1$b )

# compare estimated parameter
m2 <- cbind( sirt::mirt.wrapper.coef( mod1 )$coef, delta, tau )[,-1]
colnames(m2) <- c(  "delta.mirt", "tau.mirt", "thdim","delta.true","tau.true"   )
m2 <- cbind(m1,m2)
round( m2, 3 )
}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Beta item response model}
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
