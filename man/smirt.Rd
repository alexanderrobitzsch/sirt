%% File Name: smirt.Rd
%% File Version: 2.411

\name{smirt}
\alias{smirt}
\alias{summary.smirt}
\alias{logLik.smirt}
\alias{anova.smirt}
\alias{IRT.irfprob.smirt}
\alias{IRT.likelihood.smirt}
\alias{IRT.posterior.smirt}
\alias{IRT.modelfit.smirt}
\alias{summary.IRT.modelfit.smirt}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multidimensional Noncompensatory, Compensatory and Partially
Compensatory Item Response Model
}
\description{
This function estimates the noncompensatory and compensatory multidimensional
item response model (Bolt & Lall, 2003; Reckase, 2009) as well as
the partially compensatory item response model (Spray et al., 1990)
for dichotomous data.
}
\usage{
smirt(dat, Qmatrix, irtmodel="noncomp", est.b=NULL, est.a=NULL,
     est.c=NULL, est.d=NULL, est.mu.i=NULL, b.init=NULL, a.init=NULL,
     c.init=NULL, d.init=NULL, mu.i.init=NULL, Sigma.init=NULL,
     b.lower=-Inf, b.upper=Inf, a.lower=-Inf, a.upper=Inf,
     c.lower=-Inf, c.upper=Inf, d.lower=-Inf, d.upper=Inf,
     theta.k=seq(-6,6,len=20), theta.kDES=NULL,
     qmcnodes=0, mu.fixed=NULL, variance.fixed=NULL,  est.corr=FALSE,
     max.increment=1, increment.factor=1, numdiff.parm=0.0001,
     maxdevchange=0.1, globconv=0.001, maxiter=1000, msteps=4,
     mstepconv=0.001)

\method{summary}{smirt}(object,...)

\method{anova}{smirt}(object,...)

\method{logLik}{smirt}(object,...)

\method{IRT.irfprob}{smirt}(object,...)

\method{IRT.likelihood}{smirt}(object,...)

\method{IRT.posterior}{smirt}(object,...)

\method{IRT.modelfit}{smirt}(object,...)

\method{summary}{IRT.modelfit.smirt}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame with dichotomous item responses
}
  \item{Qmatrix}{
The Q-matrix which specifies the loadings to be estimated
}
  \item{irtmodel}{
The item response model. Options are the noncompensatory model (\code{"noncomp"}),
the compensatory model (\code{"comp"}) and
the partially compensatory model (\code{"partcomp"}).
See Details for more explanations.
}
  \item{est.b}{
An integer matrix (if \code{irtmodel="noncomp"}) or integer vector
(if \code{irtmodel="comp"}) for \eqn{b} parameters to be estimated
}
  \item{est.a}{
An integer matrix for \eqn{a} parameters to be estimated.
If \code{est.a="2PL"}, then all item loadings will be estimated and the
variances are set to one (and therefore \code{est.corr=TRUE}).
}
  \item{est.c}{
An integer vector for \eqn{c} parameters to be estimated
}
  \item{est.d}{
An integer vector for \eqn{d} parameters to be estimated
}
  \item{est.mu.i}{
An integer vector for \eqn{\mu_i} parameters to be estimated
}
  \item{b.init}{
Initial \eqn{b} coefficients. For \code{irtmodel="noncomp"} it must be a
matrix, for \code{irtmodel="comp"} it is a vector.
}
  \item{a.init}{
Initial \eqn{a} coefficients arranged in a matrix
}
  \item{c.init}{
Initial \eqn{c} coefficients
}
  \item{d.init}{
Initial \eqn{d} coefficients
}
\item{mu.i.init}{
Initial \eqn{d} coefficients
}
\item{Sigma.init}{
Initial covariance matrix \eqn{\Sigma}
}
\item{b.lower}{Lower bound for \eqn{b} parameter}
\item{b.upper}{Upper bound for \eqn{b} parameter}
\item{a.lower}{Lower bound for \eqn{a} parameter}
\item{a.upper}{Upper bound for \eqn{a} parameter}
\item{c.lower}{Lower bound for \eqn{c} parameter}
\item{c.upper}{Upper bound for \eqn{c} parameter}
\item{d.lower}{Lower bound for \eqn{d} parameter}
\item{d.upper}{Upper bound for \eqn{d} parameter}
  \item{theta.k}{
Vector of discretized trait distribution. This vector is expanded in all
dimensions by using the \code{\link[base:expand.grid]{base::expand.grid}}
function. If a user specifies a design matrix \code{theta.kDES} of transformed
\eqn{\bold{\theta}_p}
values (see Details and Examples), then \code{theta.k} must be a matrix, too.
}
  \item{theta.kDES}{
An optional design matrix. This matrix will differ from the ordinary
theta grid in case of nonlinear item response models.
}
  \item{qmcnodes}{
Number of integration nodes for quasi Monte Carlo integration (see Pan &
Thompson, 2007; Gonzales et al., 2006). Integration points are obtained by using
the function \code{\link{qmc.nodes}}. Note that when using quasi Monte Carlo nodes,
no theta design matrix \code{theta.kDES} can be specified.
See Example 1, Model 11.
}
  \item{mu.fixed}{
Matrix with fixed entries in the mean vector. By default, all means
are set to zero.
}
  \item{variance.fixed}{
Matrix (with rows and three columns) with fixed entries in the covariance matrix
(see Examples). The entry \eqn{c_{kd}} of the covariance between dimensions
  \eqn{k} and \eqn{d} is set to \eqn{c_0} iff \code{variance.fixed} has a row with
  a \eqn{k} in the first column, a \eqn{d} in the second column and the value
  \eqn{c_0} in the third column.
}
  \item{est.corr}{
Should only a correlation matrix instead of a covariance matrix
be estimated?
}
  \item{max.increment}{
Maximum increment
}
  \item{increment.factor}{
A value (larger than one) which defines the extent of the decrease of the maximum
increment of item parameters in every iteration. The maximum increment in iteration
\code{iter} is defined as \code{max.increment*increment.factor^(-iter)}
where \code{max.increment=1}. Using a value larger than 1 helps
to reach convergence in some non-converging analyses (use values of 1.01, 1.02
or even 1.05). See also Example 1 Model 2a.
}
  \item{numdiff.parm}{
Numerical differentiation parameter
}
  \item{maxdevchange}{
Convergence criterion for change in relative deviance
}
  \item{globconv}{
Global convergence criterion for parameter change
}
  \item{maxiter}{
Maximum number of iterations
}
  \item{msteps}{
Number of iterations within a M step
}
  \item{mstepconv}{
Convergence criterion within a M step
}
  \item{object}{
Object of class \code{smirt}
}
  \item{\dots}{
Further arguments to be passed
}
}
\details{
The noncompensatory item response model (\code{irtmodel="noncomp"};
e.g. Bolt & Lall, 2003) is defined as
\deqn{P(X_{pi}=1 | \bold{\theta}_p )=c_i + (d_i - c_i )
\prod_l invlogit( a_{il} q_{il}  \theta_{pl} - b_{il} ) }
where \eqn{i}, \eqn{p}, \eqn{l} denote items, persons and dimensions
respectively.

The compensatory item response model (\code{irtmodel="comp"}) is defined by
\deqn{P(X_{pi}=1 | \bold{\theta}_p )=c_i + (d_i - c_i )
invlogit( \sum_l  a_{il}  q_{il} \theta_{pl} - b_{i} ) }
Using a design matrix \code{theta.kDES} the model can be made even more general
in a model which is linear in item parameters
\deqn{P(X_{pi}=1 | \bold{\theta}_p )=c_i + (d_i - c_i )
invlogit( \sum_l  a_{il}  q_{il} t_l ( \bold{ \theta_{p} } )   - b_{i} ) }
with known functions \eqn{t_l} of the trait vector \eqn{\bold{\theta}_p}.
Fixed values of the functions \eqn{t_l} are specified in the
\eqn{\bold{\theta}_p} design matrix \code{theta.kDES}.

The partially compensatory item response model (\code{irtmodel="partcomp"})
is defined by
\deqn{P(X_{pi}=1 | \bold{\theta}_p )=c_i + (d_i - c_i )
\frac{ \exp \left( \sum_l ( a_{il}  q_{il} \theta_{pl} - b_{il} ) \right) }
{  \mu_i \prod_l ( 1 + \exp ( a_{il}  q_{il} \theta_{pl} - b_{il} ) ) +
    ( 1- \mu_i) ( 1 + \exp \left( \sum_l  ( a_{il}  q_{il} \theta_{pl} - b_{il} ) \right) )
        } }
with item parameters \eqn{\mu_i} indicating the degree of compensatory.
\eqn{\mu_i=1} indicates a noncompensatory model while \eqn{\mu_i=0}
indicates a (fully) compensatory model.

The models are estimated by an EM algorithm employing marginal maximum
likelihood.
}
\value{
A list with following entries:
\item{deviance}{Deviance}
\item{ic}{Information criteria}
\item{item}{Data frame with item parameters}
\item{person}{Data frame with person parameters. It includes
    the person mean of all item responses (\code{M}; percentage correct
    of all non-missing items), the
EAP estimate and its corresponding standard error for all dimensions
(\code{EAP} and \code{SE.EAP}) and the maximum likelihood estimate
as well as the mode of the posterior distribution (\code{MLE} and \code{MAP}).}
\item{EAP.rel}{EAP reliability}
\item{mean.trait}{Means of trait}
\item{sd.trait}{Standard deviations of trait}
\item{Sigma}{Trait covariance matrix}
\item{cor.trait}{Trait correlation matrix}
\item{b}{Matrix (vector) of \eqn{b} parameters}
\item{se.b}{Matrix (vector) of standard errors \eqn{b} parameters}
\item{a}{Matrix of \eqn{a} parameters}
\item{se.a}{Matrix of standard errors of \eqn{a} parameters}
\item{c}{Vector of \eqn{c} parameters}
\item{se.c}{Vector of standard errors of \eqn{c} parameters}
\item{d}{Vector of \eqn{d} parameters}
\item{se.d}{Vector of standard errors of \eqn{d} parameters}
\item{mu.i}{Vector of \eqn{\mu_i} parameters}
\item{se.mu.i}{Vector of standard errors of \eqn{\mu_i} parameters}
\item{f.yi.qk}{Individual likelihood}
\item{f.qk.yi}{Individual posterior}
\item{probs}{Probabilities of item response functions evaluated at \code{theta.k}}
\item{n.ik}{Expected counts}
\item{iter}{Number of iterations}
\item{dat2}{Processed data set}
\item{dat2.resp}{Data set of response indicators}
\item{I}{Number of items}
\item{D}{Number of dimensions}
\item{K}{Maximum item response score}
\item{theta.k}{Used theta integration grid}
\item{pi.k}{Distribution function evaluated at \code{theta.k}}
\item{irtmodel}{Used IRT model}
\item{Qmatrix}{Used Q-matrix}
}
\references{
Bolt, D. M., & Lall, V. F. (2003). Estimation of compensatory and
noncompensatory multidimensional item response models using Markov chain
Monte Carlo. \emph{Applied Psychological Measurement, 27}, 395-414.

Gonzalez, J., Tuerlinckx, F., De Boeck, P., & Cools, R. (2006).
Numerical integration in logistic-normal models.
\emph{Computational Statistics & Data Analysis, 51}, 1535-1548.

Pan, J., & Thompson, R. (2007). Quasi-Monte Carlo estimation in
generalized linear mixed models. \emph{Computational Statistics &
Data Analysis, 51}, 5765-5775.

Reckase, M. (2009). \emph{Multidimensional item response theory}.
New York: Springer. \doi{10.1007/978-0-387-89976-3}

Spray, J. A., Davey, T. C., Reckase, M. D., Ackerman, T. A., & Carlson, J. E. (1990).
\emph{Comparison of two logistic multidimensional item response theory models}.
ACT Research Report No. ACT-RR-ONR-90-8.
}


%\author{
%Alexander Robitzsch
%}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See the \code{\link[mirt:mirt]{mirt::mirt}} and
\code{itemtype="partcomp"} for estimating noncompensatory item response models
using the \pkg{mirt} package. See also \code{\link[mirt:mixedmirt]{mirt::mixedmirt}}.

Other multidimensional IRT models can also be estimated
with \code{\link{rasch.mml2}} and \code{\link{rasch.mirtlc}}.

See \code{\link[CDM:itemfit.sx2]{itemfit.sx2}} (\pkg{CDM}) for item fit
statistics.

See also the \pkg{mirt} and \pkg{TAM} packages for estimation of
compensatory multidimensional item response models.
}
\examples{
#############################################################################
## EXAMPLE 1: Noncompensatory and compensatory IRT models
#############################################################################
set.seed(997)

# (1) simulate data from a two-dimensional noncompensatory
#     item response model
#   -> increase number of iterations in all models!

N <- 1000    # number of persons
I <- 10        # number of items
theta0 <- rnorm( N, sd=1 )
theta1 <- theta0 + rnorm(N, sd=.7 )
theta2 <- theta0 + rnorm(N, sd=.7 )
Q <- matrix( 1, nrow=I,ncol=2 )
Q[ 1:(I/2), 2 ] <- 0
Q[ I,1] <- 0
b <- matrix( rnorm( I*2 ), I, 2 )
a <- matrix( 1, I, 2 )

# simulate data
prob <- dat <- matrix(0, nrow=N, ncol=I )
for (ii in 1:I){
prob[,ii] <- ( stats::plogis( theta1 - b[ii,1]  ) )^Q[ii,1]
prob[,ii] <- prob[,ii] * ( stats::plogis( theta2 - b[ii,2]  ) )^Q[ii,2]
            }
dat[ prob > matrix( stats::runif( N*I),N,I) ] <- 1
colnames(dat) <- paste0("I",1:I)

#***
# Model 1: Noncompensatory 1PL model
mod1 <- sirt::smirt(dat, Qmatrix=Q, maxiter=10 ) # change number of iterations
summary(mod1)

\dontrun{
#***
# Model 2: Noncompensatory 2PL model
mod2 <- sirt::smirt(dat,Qmatrix=Q, est.a="2PL", maxiter=15 )
summary(mod2)

# Model 2a: avoid convergence problems with increment.factor
mod2a <- sirt::smirt(dat,Qmatrix=Q, est.a="2PL", maxiter=30, increment.factor=1.03)
summary(mod2a)

#***
# Model 3: some fixed c and d parameters different from zero or one
c.init <- rep(0,I)
c.init[ c(3,7)] <- .2
d.init <- rep(1,I)
d.init[c(4,8)] <- .95
mod3 <- sirt::smirt( dat, Qmatrix=Q, c.init=c.init, d.init=d.init )
summary(mod3)

#***
# Model 4: some estimated c and d parameters (in parameter groups)
est.c <- c.init <- rep(0,I)
c.estpars <- c(3,6,7)
c.init[ c.estpars ] <- .2
est.c[c.estpars] <- 1
est.d <- rep(0,I)
d.init <- rep(1,I)
d.estpars <- c(6,9)
d.init[ d.estpars ] <- .95
est.d[ d.estpars ] <- d.estpars   # different d parameters
mod4 <- sirt::smirt(dat,Qmatrix=Q, est.c=est.c, c.init=c.init,
            est.d=est.d, d.init=d.init  )
summary(mod4)

#***
# Model 5: Unidimensional 1PL model
Qmatrix <- matrix( 1, nrow=I, ncol=1 )
mod5 <- sirt::smirt( dat, Qmatrix=Qmatrix )
summary(mod5)

#***
# Model 6: Unidimensional 2PL model
mod6 <- sirt::smirt( dat, Qmatrix=Qmatrix, est.a="2PL" )
summary(mod6)

#***
# Model 7: Compensatory model with between item dimensionality
# Note that the data is simulated under the noncompensatory condition
# Therefore Model 7 should have a worse model fit than Model 1
Q1 <- Q
Q1[ 6:10, 1] <- 0
mod7 <- sirt::smirt(dat,Qmatrix=Q1, irtmodel="comp", maxiter=30)
summary(mod7)

#***
# Model 8: Compensatory model with within item dimensionality
#         assuming zero correlation between dimensions
variance.fixed <- as.matrix( cbind( 1,2,0) )
# set the covariance between the first and second dimension to zero
mod8 <- sirt::smirt(dat,Qmatrix=Q, irtmodel="comp", variance.fixed=variance.fixed,
            maxiter=30)
summary(mod8)

#***
# Model 8b: 2PL model with starting values for a and b parameters
b.init <- rep(0,10)  # set all item difficulties initially to zero
# b.init <- NULL
a.init <- Q       # initialize a.init with Q-matrix
# provide starting values for slopes of first three items on Dimension 1
a.init[1:3,1] <- c( .55, .32, 1.3)

mod8b <- sirt::smirt(dat,Qmatrix=Q, irtmodel="comp", variance.fixed=variance.fixed,
              b.init=b.init, a.init=a.init, maxiter=20, est.a="2PL" )
summary(mod8b)

#***
# Model 9: Unidimensional model with quadratic item response functions
# define theta
theta.k <- seq( - 6, 6, len=15 )
theta.k <- as.matrix( theta.k, ncol=1 )
# define design matrix
theta.kDES <- cbind( theta.k[,1], theta.k[,1]^2 )
# define Q-matrix
Qmatrix <- matrix( 0, I, 2 )
Qmatrix[,1] <- 1
Qmatrix[ c(3,6,7), 2 ] <- 1
colnames(Qmatrix) <- c("F1", "F1sq" )
# estimate model
mod9 <- sirt::smirt(dat,Qmatrix=Qmatrix, maxiter=50, irtmodel="comp",
           theta.k=theta.k, theta.kDES=theta.kDES, est.a="2PL" )
summary(mod9)

#***
# Model 10: Two-dimensional item response model with latent interaction
#           between dimensions
theta.k <- seq( - 6, 6, len=15 )
theta.k <- expand.grid( theta.k, theta.k )    # expand theta to 2 dimensions
# define design matrix
theta.kDES <- cbind( theta.k, theta.k[,1]*theta.k[,2] )
# define Q-matrix
Qmatrix <- matrix( 0, I, 3 )
Qmatrix[,1] <- 1
Qmatrix[ 6:10, c(2,3) ] <- 1
colnames(Qmatrix) <- c("F1", "F2", "F1iF2" )
# estimate model
mod10 <- sirt::smirt(dat,Qmatrix=Qmatrix,irtmodel="comp", theta.k=theta.k,
            theta.kDES=theta.kDES, est.a="2PL" )
summary(mod10)

#****
# Model 11: Example Quasi Monte Carlo integration
Qmatrix <- matrix( 1, I, 1 )
mod11 <- sirt::smirt( dat, irtmodel="comp", Qmatrix=Qmatrix, qmcnodes=1000 )
summary(mod11)

#############################################################################
## EXAMPLE 2: Dataset Reading data.read
##            Multidimensional models for dichotomous data
#############################################################################

data(data.read)
dat <- data.read
I <- ncol(dat)    # number of items

#***
# Model 1: 3-dimensional 2PL model

# define Q-matrix
Qmatrix <- matrix(0,nrow=I,ncol=3)
Qmatrix[1:4,1] <- 1
Qmatrix[5:8,2] <- 1
Qmatrix[9:12,3] <- 1

# estimate model
mod1 <- sirt::smirt( dat, Qmatrix=Qmatrix, irtmodel="comp", est.a="2PL",
            qmcnodes=1000, maxiter=20)
summary(mod1)

#***
# Model 2: 3-dimensional Rasch model
mod2 <- sirt::smirt( dat, Qmatrix=Qmatrix, irtmodel="comp",
              qmcnodes=1000, maxiter=20)
summary(mod2)

#***
# Model 3: 3-dimensional 2PL model with uncorrelated dimensions
# fix entries in variance matrix
variance.fixed <- cbind( c(1,1,2), c(2,3,3), 0 )
# set the following covariances to zero: cov[1,2]=cov[1,3]=cov[2,3]=0

# estimate model
mod3 <- sirt::smirt( dat, Qmatrix=Qmatrix, irtmodel="comp", est.a="2PL",
             variance.fixed=variance.fixed, qmcnodes=1000, maxiter=20)
summary(mod3)

#***
# Model 4: Bifactor model with one general factor (g) and
#          uncorrelated specific factors

# define a new Q-matrix
Qmatrix1 <- cbind( 1, Qmatrix )
# uncorrelated factors
variance.fixed <- cbind( c(1,1,1,2,2,3), c(2,3,4,3,4,4), 0 )
# The first dimension refers to the general factors while the other
# dimensions refer to the specific factors.
# The specification means that:
# Cov[1,2]=Cov[1,3]=Cov[1,4]=Cov[2,3]=Cov[2,4]=Cov[3,4]=0

# estimate model
mod4 <- sirt::smirt( dat, Qmatrix=Qmatrix1, irtmodel="comp", est.a="2PL",
             variance.fixed=variance.fixed, qmcnodes=1000, maxiter=20)
summary(mod4)

#############################################################################
## EXAMPLE 3: Partially compensatory model
#############################################################################

#**** simulate data
set.seed(7656)
I <- 10         # number of items
N <- 2000        # number of subjects
Q <- matrix( 0, 3*I,2)  # Q-matrix
Q[1:I,1] <- 1
Q[1:I + I,2] <- 1
Q[1:I + 2*I,1:2] <- 1
b <- matrix( stats::runif( 3*I *2, -2, 2 ), nrow=3*I, 2 )
b <- b*Q
b <- round( b, 2 )
mui <- rep(0,3*I)
mui[ seq(2*I+1, 3*I) ] <- 0.65
# generate data
dat <- matrix( NA, N, 3*I )
colnames(dat) <- paste0("It", 1:(3*I) )
# simulate item responses
library(mvtnorm)
theta <- mvtnorm::rmvnorm(N, mean=c(0,0), sigma=matrix( c( 1.2, .6,.6,1.6),2, 2 ) )
for (ii in 1:(3*I)){
    # define probability
    tmp1 <- exp( theta[,1] * Q[ii,1] - b[ii,1] +  theta[,2] * Q[ii,2] - b[ii,2] )
    # non-compensatory model
    nco1 <- ( 1 + exp( theta[,1] * Q[ii,1] - b[ii,1] ) ) *
                  ( 1 + exp( theta[,2] * Q[ii,2] - b[ii,2] ) )
    co1 <- ( 1 + tmp1 )
    p1 <- tmp1 / ( mui[ii] * nco1 + ( 1 - mui[ii] )*co1 )
    dat[,ii] <- 1 * ( stats::runif(N) < p1 )
}

#*** Model 1: Joint mu.i parameter for all items
est.mu.i <- rep(0,3*I)
est.mu.i[ seq(2*I+1,3*I)] <- 1
mod1 <- sirt::smirt( dat, Qmatrix=Q, irtmodel="partcomp", est.mu.i=est.mu.i)
summary(mod1)

#*** Model 2: Separate mu.i parameter for all items
est.mu.i[ seq(2*I+1,3*I)] <- 1:I
mod2 <- sirt::smirt( dat, Qmatrix=Q, irtmodel="partcomp", est.mu.i=est.mu.i)
summary(mod2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Multidimensional item response model}
%% \keyword{summary}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
