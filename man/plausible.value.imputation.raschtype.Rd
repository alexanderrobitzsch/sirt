%% File Name: plausible.value.imputation.raschtype.Rd
%% File Version: 1.34

\name{plausible.value.imputation.raschtype}
\alias{plausible.value.imputation.raschtype}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plausible Value Imputation in Generalized Logistic Item
Response Model
}
\description{
This function performs unidimensional plausible value imputation
(Adams & Wu, 2007; Mislevy, 1991).
}
\usage{
plausible.value.imputation.raschtype(data=NULL, f.yi.qk=NULL, X,
   Z=NULL, beta0=rep(0, ncol(X)), sig0=1, b=rep(1, ncol(X)),
   a=rep(1, length(b)), c=rep(0, length(b)), d=1+0*b,
   alpha1=0, alpha2=0, theta.list=seq(-5, 5, len=50),
   cluster=NULL, iter, burnin, nplausible=1, printprogress=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
An \eqn{N \times I} data frame of dichotomous responses
}
\item{f.yi.qk}{
An optional matrix which contains the individual likelihood.
  This matrix is produced by \code{\link{rasch.mml2}} or
  \code{\link{rasch.copula2}}. The use of this argument allows the
  estimation of the latent regression model independent of the
  parameters of the used item response model.
    }
\item{X}{
A matrix of individual covariates for the latent
regression of \eqn{\theta} on \eqn{X}
}
\item{Z}{
A matrix of individual covariates for the regression
of individual residual variances on \eqn{Z}
}
  \item{beta0}{
Initial vector of regression coefficients
}
  \item{sig0}{
Initial vector of coefficients for the
variance heterogeneity model
}
  \item{b}{
Vector of item difficulties. It must not be provided
if the individual likelihood \code{f.yi.qk} is
specified.
}
  \item{a}{
Optional vector of item slopes
}
  \item{c}{
Optional vector of lower item asymptotes
}
  \item{d}{
Optional vector of upper item asymptotes
}
  \item{alpha1}{
Parameter \eqn{\alpha_1} in generalized
item response model
}
  \item{alpha2}{
Parameter \eqn{\alpha_2} in generalized
item response model
}
  \item{theta.list}{
Vector of theta values at which the ability
distribution should be evaluated
}
  \item{cluster}{
Cluster identifier (e.g. schools or classes)
for including theta means in the plausible imputation.
}
  \item{iter}{
Number of iterations
}
\item{burnin}{
Number of burn-in iterations for plausible value imputation
}
  \item{nplausible}{
Number of plausible values
}
  \item{printprogress}{
A logical indicated whether iteration progress should be displayed
at the console.
}
}

\details{
Plausible values are drawn from the latent regression model with
heterogeneous variances:
\deqn{\theta_p=X_p \beta + \epsilon_p  \quad, \quad
\epsilon_p \sim N( 0, \sigma_p^2 ) \quad, \quad
\log( \sigma_p )=Z_p \gamma + \nu_p }
}


\value{
A list with following entries:
\item{coefs.X}{Sampled regression coefficients for covariates \eqn{X}}
\item{coefs.Z}{Sampled coefficients for modeling variance heterogeneity
for covariates \eqn{Z}}
\item{pvdraws}{Matrix with drawn plausible values}
\item{posterior}{Posterior distribution from last iteration}
\item{EAP}{Individual EAP estimate}
\item{SE.EAP}{Standard error of the EAP estimate}
\item{pv.indexes}{Index of iterations for which plausible
values were drawn}
}

\references{
Adams, R., & Wu. M. (2007). The mixed-coefficients multinomial
logit model: A generalized form of the Rasch model.
In M. von Davier & C. H. Carstensen: \emph{Multivariate and Mixture
Distribution Rasch Models: Extensions and Applications} (pp. 57-76).
New York: Springer.

Mislevy, R. J. (1991). Randomization-based inference about latent
variables from complex samples. \emph{Psychometrika, 56}, 177-196.
}


%\author{
%Alexander Robitzsch
%}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For estimating the latent regression model see
\code{\link{latent.regression.em.raschtype}}.
}



\examples{
#############################################################################
# EXAMPLE 1: Rasch model with covariates
#############################################################################

set.seed(899)
I <- 21     # number of items
b <- seq(-2,2, len=I)   # item difficulties
n <- 2000       # number of students

# simulate theta and covariates
theta <- stats::rnorm( n )
x <- .7 * theta + stats::rnorm( n, .5 )
y <- .2 * x+ .3*theta + stats::rnorm( n, .4 )
dfr <- data.frame( theta, 1, x, y )

# simulate Rasch model
dat1 <- sirt::sim.raschtype( theta=theta, b=b )

# Plausible value draws
pv1 <- sirt::plausible.value.imputation.raschtype(data=dat1, X=dfr[,-1], b=b,
            nplausible=3, iter=10, burnin=5)
# estimate linear regression based on first plausible value
mod1 <- stats::lm( pv1$pvdraws[,1] ~ x+y )
summary(mod1)
  ##               Estimate Std. Error t value Pr(>|t|)
  ##   (Intercept) -0.27755    0.02121  -13.09   <2e-16 ***
  ##   x            0.40483    0.01640   24.69   <2e-16 ***
  ##   y            0.20307    0.01822   11.15   <2e-16 ***

# true regression estimate
summary( stats::lm( theta ~ x + y ) )
  ## Coefficients:
  ##             Estimate Std. Error t value Pr(>|t|)
  ## (Intercept) -0.27821    0.01984  -14.02   <2e-16 ***
  ## x            0.40747    0.01534   26.56   <2e-16 ***
  ## y            0.18189    0.01704   10.67   <2e-16 ***

\dontrun{
#############################################################################
# EXAMPLE 2: Classical test theory, homogeneous regression variance
#############################################################################

set.seed(899)
n <- 3000       # number of students
x <- round( stats::runif( n, 0,1 ) )
y <- stats::rnorm(n)
# simulate true score theta
theta <- .4*x + .5 * y + stats::rnorm(n)
# simulate observed score by adding measurement error
sig.e <- rep( sqrt(.40), n )
theta_obs <- theta + stats::rnorm( n, sd=sig.e)

# define theta grid for evaluation of density
theta.list <- mean(theta_obs) + stats::sd(theta_obs) * seq( - 5, 5, length=21)
# compute individual likelihood
f.yi.qk <- stats::dnorm( outer( theta_obs, theta.list, "-" ) / sig.e )
f.yi.qk <- f.yi.qk / rowSums(f.yi.qk)
# define covariates
X <- cbind( 1, x, y )
# draw plausible values
mod2 <- sirt::plausible.value.imputation.raschtype( f.yi.qk=f.yi.qk,
                  theta.list=theta.list, X=X, iter=10, burnin=5)

# linear regression
mod1 <- stats::lm( mod2$pvdraws[,1] ~ x+y )
summary(mod1)
  ##             Estimate Std. Error t value Pr(>|t|)
  ## (Intercept) -0.01393    0.02655  -0.525      0.6
  ## x            0.35686    0.03739   9.544   <2e-16 ***
  ## y            0.53759    0.01872  28.718   <2e-16 ***

# true regression model
summary( stats::lm( theta ~ x + y ) )
  ##             Estimate Std. Error t value Pr(>|t|)
  ## (Intercept) 0.002931   0.026171   0.112    0.911
  ## x           0.359954   0.036864   9.764   <2e-16 ***
  ## y           0.509073   0.018456  27.584   <2e-16 ***

#############################################################################
# EXAMPLE 3: Classical test theory, heterogeneous regression variance
#############################################################################

set.seed(899)
n <- 5000       # number of students
x <- round( stats::runif( n, 0,1 ) )
y <- stats::rnorm(n)
# simulate true score theta
theta <- .4*x + .5 * y + stats::rnorm(n) * ( 1 - .4 * x )
# simulate observed score by adding measurement error
sig.e <- rep( sqrt(.40), n )
theta_obs <- theta + stats::rnorm( n, sd=sig.e)

# define theta grid for evaluation of density
theta.list <- mean(theta_obs) + stats::sd(theta_obs) * seq( - 5, 5, length=21)
# compute individual likelihood
f.yi.qk <- stats::dnorm( outer( theta_obs, theta.list, "-" ) / sig.e )
f.yi.qk <- f.yi.qk / rowSums(f.yi.qk)
# define covariates
X <- cbind( 1, x, y )
# draw plausible values (assuming variance homogeneity)
mod3a <- sirt::plausible.value.imputation.raschtype( f.yi.qk=f.yi.qk,
                  theta.list=theta.list, X=X, iter=10, burnin=5)
# draw plausible values (assuming variance heterogeneity)
#  -> include predictor Z
mod3b <- sirt::plausible.value.imputation.raschtype( f.yi.qk=f.yi.qk,
                  theta.list=theta.list, X=X, Z=X, iter=10, burnin=5)

# investigate variance of theta conditional on x
res3 <- sapply( 0:1, FUN=function(vv){
        c( stats::var(theta[x==vv]), stats::var(mod3b$pvdraw[x==vv,1]),
              stats::var(mod3a$pvdraw[x==vv,1]))})
rownames(res3) <- c("true", "pv(hetero)", "pv(homog)" )
colnames(res3) <- c("x=0","x=1")
  ## > round( res3, 2 )
  ##             x=0  x=1
  ## true       1.30 0.58
  ## pv(hetero) 1.29 0.55
  ## pv(homog)  1.06 0.77
## -> assuming heteroscedastic variances recovers true conditional variance
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Plausible values}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
