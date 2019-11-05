%% File Name: latent.regression.em.raschtype.Rd
%% File Version: 1.481

\name{latent.regression.em.raschtype}
\alias{latent.regression.em.raschtype}
\alias{latent.regression.em.normal}
\alias{summary.latent.regression}


\title{
Latent Regression Model for the Generalized
Logistic Item Response Model and the Linear Model for Normal Responses
}
\description{
This function estimates a unidimensional latent regression model if
a likelihood is specified, parameters from the generalized
item response model (Stukel, 1988) or a mean and a standard error
estimate for individual scores is provided as input.
Item parameters are treated as fixed in the estimation.
}

\usage{
latent.regression.em.raschtype(data=NULL, f.yi.qk=NULL, X,
    weights=rep(1, nrow(X)), beta.init=rep(0,ncol(X)),
    sigma.init=1, b=rep(0,ncol(X)), a=rep(1,length(b)),
    c=rep(0, length(b)), d=rep(1, length(b)), alpha1=0, alpha2=0,
    max.parchange=1e-04, theta.list=seq(-5, 5, len=20),
    maxiter=300, progress=TRUE )

latent.regression.em.normal(y, X, sig.e, weights=rep(1, nrow(X)),
    beta.init=rep(0, ncol(X)), sigma.init=1, max.parchange=1e-04,
    maxiter=300, progress=TRUE)

\method{summary}{latent.regression}(object,...)
}

\arguments{
  \item{data}{ An \eqn{N \times I} data frame of dichotomous item responses.
  If no data frame is supplied, then a user can input the individual
  likelihood \code{f.yi.qk}.
   }
  \item{f.yi.qk}{An optional matrix which contains the individual likelihood.
  This matrix is produced by \code{\link{rasch.mml2}} or
  \code{\link{rasch.copula2}}. The use of this argument allows the
  estimation of the latent regression model independent of the
  parameters of the used item response model.
    }
  \item{X}{An \eqn{N \times K} matrix of \eqn{K} covariates in the latent
    regression model. Note that the intercept (i.e. a vector of ones)
    must be included in \code{X}.
}
  \item{weights}{
Student weights (optional).
}
  \item{beta.init}{
Initial regression coefficients (optional).
}
\item{sigma.init}{
Initial residual standard deviation (optional).
}
  \item{b}{
Item difficulties (optional). They must only be provided
if the likelihood \code{f.yi.qk} is not given as an input.
}
  \item{a}{
Item discriminations (optional).
}
  \item{c}{
Guessing parameter (lower asymptotes) (optional).
}
  \item{d}{
One minus slipping parameter (upper asymptotes) (optional).
}
  \item{alpha1}{
Upper tail parameter \eqn{\alpha_1} in the generalized logistic item response model.
Default is 0.
}
  \item{alpha2}{
Lower tail parameter \eqn{\alpha_2} parameter in the generalized
logistic item response model. Default is 0.
}
  \item{max.parchange}{
Maximum change in regression parameters
}
  \item{theta.list}{
Grid of person ability where theta is evaluated
}
  \item{maxiter}{
Maximum number of iterations
}
\item{progress}{An optional logical indicating whether computation
progress should be displayed.}
\item{y}{Individual scores}
\item{sig.e}{Standard errors for individual scores}
  \item{object}{
Object of class \code{latent.regression}
}
  \item{\dots}{
Further arguments to be passed
}

}
\details{
In the output \emph{Regression Parameters}
the fraction of missing information (\code{fmi}) is reported
which is the increase of variance in regression
parameter estimates because ability is defined as
a latent variable. The effective sample size \code{pseudoN.latent}
corresponds to a sample size when the ability would be
available with a reliability of one.
}


\value{
A list with following entries
  \item{iterations}{Number of iterations needed}
  \item{maxiter}{Maximal number of iterations}
  \item{max.parchange}{Maximum change in parameter estimates}
  \item{coef}{Coefficients}
  \item{summary.coef}{Summary of regression coefficients}
  \item{sigma}{Estimate of residual standard deviation}
  \item{vcov.simple}{Covariance parameters of estimated parameters
      (simplified version)}
  \item{vcov.latent}{Covariance parameters of estimated parameters
                        which accounts for latent ability}
\item{post}{Individual posterior distribution}
\item{EAP}{Individual EAP estimates}
\item{SE.EAP}{Standard error estimates of EAP}
  \item{explvar}{Explained variance in latent regression}
  \item{totalvar}{Total variance in latent regression}
  \item{rsquared}{Explained variance \eqn{R^2} in latent regression}
}

\references{
Adams, R., & Wu. M. (2007). The mixed-coefficients multinomial
logit model: A generalized form of the Rasch model.
In M. von Davier & C. H. Carstensen (Eds.). \emph{Multivariate and mixture
distribution Rasch models: Extensions and applications} (pp. 57-76).
New York: Springer.
\doi{10.1007/978-0-387-49839-3_4}

Mislevy, R. J. (1991). Randomization-based inference about latent variables
from complex samples. \emph{Psychometrika, 56}(2), 177-196.
\doi{10.1007/BF02294457}

Stukel, T. A. (1988). Generalized logistic models.
\emph{Journal of the American Statistical Association, 83}(402), 426-431.
\doi{10.1080/01621459.1988.10478613}
}

%\author{
%Alexander Robitzsch
%}


\note{
Using the defaults in \code{a}, \code{c}, \code{d},
\code{alpha1} and \code{alpha2} corresponds to the
Rasch model.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{plausible.value.imputation.raschtype}}
for plausible value imputation of generalized logistic
item type models.
}

\examples{
#############################################################################
#  EXAMPLE 1: PISA Reading | Rasch model for dichotomous data
#############################################################################

data(data.pisaRead, package="sirt")
dat <- data.pisaRead$data
items <- grep("R", colnames(dat))
# define matrix of covariates
X <- cbind( 1, dat[, c("female","hisei","migra" ) ] )

#***
# Model 1: Latent regression model in the Rasch model
# estimate Rasch model
mod1 <- sirt::rasch.mml2( dat[,items] )
# latent regression model
lm1 <- sirt::latent.regression.em.raschtype( data=dat[,items ], X=X, b=mod1$item$b )

\dontrun{
#***
# Model 2: Latent regression with generalized link function
# estimate alpha parameters for link function
mod2 <- sirt::rasch.mml2( dat[,items], est.alpha=TRUE)
# use model estimated likelihood for latent regression model
lm2 <- sirt::latent.regression.em.raschtype( f.yi.qk=mod2$f.yi.qk,
            X=X, theta.list=mod2$theta.k)

#***
# Model 3: Latent regression model based on Rasch copula model
testlets <- paste( data.pisaRead$item$testlet)
itemclusters <- match( testlets, unique(testlets) )
# estimate Rasch copula model
mod3 <- sirt::rasch.copula2( dat[,items], itemcluster=itemclusters )
# use model estimated likelihood for latent regression model
lm3 <- sirt::latent.regression.em.raschtype( f.yi.qk=mod3$f.yi.qk,
                X=X, theta.list=mod3$theta.k)

#############################################################################
# EXAMPLE 2: Simulated data according to the Rasch model
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

# estimate latent regression
mod <- sirt::latent.regression.em.raschtype( data=dat1, X=dfr[,-1], b=b )
  ## Regression Parameters
  ##
  ##        est se.simple     se        t p   beta    fmi N.simple pseudoN.latent
  ## X1 -0.2554    0.0208 0.0248 -10.2853 0 0.0000 0.2972     2000       1411.322
  ## x   0.4113    0.0161 0.0193  21.3037 0 0.4956 0.3052     2000       1411.322
  ## y   0.1715    0.0179 0.0213   8.0438 0 0.1860 0.2972     2000       1411.322
  ##
  ## Residual Variance=0.685
  ## Explained Variance=0.3639
  ## Total Variance=1.049
  ##                 R2=0.3469

# compare with linear model (based on true scores)
summary( stats::lm( theta  ~ x + y, data=dfr ) )
  ## Coefficients:
  ##             Estimate Std. Error t value Pr(>|t|)
  ## (Intercept) -0.27821    0.01984  -14.02   <2e-16 ***
  ## x            0.40747    0.01534   26.56   <2e-16 ***
  ## y            0.18189    0.01704   10.67   <2e-16 ***
  ## ---
  ##
  ## Residual standard error: 0.789 on 1997 degrees of freedom
  ## Multiple R-squared: 0.3713,     Adjusted R-squared: 0.3707

#***********
# define guessing parameters (lower asymptotes) and
# upper asymptotes ( 1 minus slipping parameters)
cI <- rep(.2, I)        # all items get a guessing parameter of .2
cI[ c(7,9) ] <- .25     # 7th and 9th get a guessing parameter of .25
dI <- rep( .95, I )    # upper asymptote of .95
dI[ c(7,11) ] <- 1        # 7th and 9th item have an asymptote of 1

# latent regression model
mod1 <- sirt::latent.regression.em.raschtype( data=dat1, X=dfr[,-1],
           b=b, c=cI, d=dI    )
  ## Regression Parameters
  ##
  ##        est se.simple     se        t p   beta    fmi N.simple pseudoN.latent
  ## X1 -0.7929    0.0243 0.0315 -25.1818 0 0.0000 0.4044     2000       1247.306
  ## x   0.5025    0.0188 0.0241  20.8273 0 0.5093 0.3936     2000       1247.306
  ## y   0.2149    0.0209 0.0266   8.0850 0 0.1960 0.3831     2000       1247.306
  ##
  ## Residual Variance=0.9338
  ## Explained Variance=0.5487
  ## Total Variance=1.4825
  ##                 R2=0.3701

#############################################################################
# EXAMPLE 3: Measurement error in dependent variable
#############################################################################

set.seed(8766)
N <- 4000       # number of persons
X <- stats::rnorm(N)           # independent variable
Z <- stats::rnorm(N)           # independent variable
y <- .45 * X + .25 * Z + stats::rnorm(N)   # dependent variable true score
sig.e <- stats::runif( N, .5, .6 )       # measurement error standard deviation
yast <- y + stats::rnorm( N, sd=sig.e ) # dependent variable measured with error

#****
# Model 1: Estimation with latent.regression.em.raschtype using
#          individual likelihood
# define theta grid for evaluation of density
theta.list <- mean(yast) + stats::sd(yast) * seq( - 5, 5, length=21)
# compute individual likelihood
f.yi.qk <- stats::dnorm( outer( yast, theta.list, "-" ) / sig.e )
f.yi.qk <- f.yi.qk / rowSums(f.yi.qk)
# define predictor matrix
X1 <- as.matrix(data.frame( "intercept"=1, "X"=X, "Z"=Z ))

# latent regression model
res <- sirt::latent.regression.em.raschtype( f.yi.qk=f.yi.qk,
                    X=X1, theta.list=theta.list)
  ##   Regression Parameters
  ##
  ##                est se.simple     se       t      p   beta    fmi N.simple pseudoN.latent
  ##   intercept 0.0112    0.0157 0.0180  0.6225 0.5336 0.0000 0.2345     4000       3061.998
  ##   X         0.4275    0.0157 0.0180 23.7926 0.0000 0.3868 0.2350     4000       3061.998
  ##   Z         0.2314    0.0156 0.0178 12.9868 0.0000 0.2111 0.2349     4000       3061.998
  ##
  ##   Residual Variance=0.9877
  ##   Explained Variance=0.2343
  ##   Total Variance=1.222
  ##                   R2=0.1917

#****
# Model 2: Estimation with latent.regression.em.normal
res2 <- sirt::latent.regression.em.normal( y=yast, sig.e=sig.e, X=X1)
  ##   Regression Parameters
  ##
  ##                est se.simple     se       t      p   beta    fmi N.simple pseudoN.latent
  ##   intercept 0.0112    0.0157 0.0180  0.6225 0.5336 0.0000 0.2345     4000       3062.041
  ##   X         0.4275    0.0157 0.0180 23.7927 0.0000 0.3868 0.2350     4000       3062.041
  ##   Z         0.2314    0.0156 0.0178 12.9870 0.0000 0.2111 0.2349     4000       3062.041
  ##
  ##   Residual Variance=0.9877
  ##   Explained Variance=0.2343
  ##   Total Variance=1.222
  ##                   R2=0.1917

  ## -> Results between Model 1 and Model 2 are identical because they use
  ##    the same input.

#***
# Model 3: Regression model based on true scores y
mod3 <- stats::lm( y ~ X + Z )
summary(mod3)
  ##   Coefficients:
  ##               Estimate Std. Error t value Pr(>|t|)
  ##   (Intercept)  0.02364    0.01569   1.506    0.132
  ##   X            0.42401    0.01570  27.016   <2e-16 ***
  ##   Z            0.23804    0.01556  15.294   <2e-16 ***
  ##   Residual standard error: 0.9925 on 3997 degrees of freedom
  ##   Multiple R-squared:  0.1923,    Adjusted R-squared:  0.1919
  ##   F-statistic: 475.9 on 2 and 3997 DF,  p-value: < 2.2e-16

#***
# Model 4: Regression model based on observed scores yast
mod4 <- stats::lm( yast ~ X + Z )
summary(mod4)
  ##   Coefficients:
  ##               Estimate Std. Error t value Pr(>|t|)
  ##   (Intercept)  0.01101    0.01797   0.613     0.54
  ##   X            0.42716    0.01797  23.764   <2e-16 ***
  ##   Z            0.23174    0.01783  13.001   <2e-16 ***
  ##   Residual standard error: 1.137 on 3997 degrees of freedom
  ##   Multiple R-squared:  0.1535,    Adjusted R-squared:  0.1531
  ##   F-statistic: 362.4 on 2 and 3997 DF,  p-value: < 2.2e-16
}
}
%% \keyword{Latent regression model}
%% \keyword{Generalized logistic item response model}
%% \keyword{Plausible values}
%% \keyword{summary}
