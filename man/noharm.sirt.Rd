%% File Name: noharm.sirt.Rd
%% File Version: 0.41

\name{noharm.sirt}
\alias{noharm.sirt}
\alias{summary.noharm.sirt}
% \alias{plot.noharm.sirt}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
NOHARM Model in \R
}

\description{
The function is an \R implementation of the normal ogive harmonic analysis
robust method (the NOHARM model; McDonald, 1997). Exploratory and confirmatory
multidimensional item response models for dichotomous data using the
probit link function can be estimated. Lower asymptotes (guessing parameters)
and upper asymptotes (One minus slipping parameters) can be provided as
fixed values.
}

\usage{
noharm.sirt(dat, weights=NULL, Fval=NULL, Fpatt=NULL, Pval=NULL,
   Ppatt=NULL, Psival=NULL, Psipatt=NULL, dimensions=NULL,
   lower=rep(0, ncol(dat)), upper=rep(1, ncol(dat)), wgtm=NULL, modesttype=1,
   pos.loading=FALSE, pos.variance=FALSE,  pos.residcorr=FALSE,
   maxiter=1000, conv=10^(-6), increment.factor=1.01, reliability=TRUE)

\method{summary}{noharm.sirt}(object, file=NULL, \dots)
}

%% \method{plot}{noharm.sirt}(x, what="est", efa.load.min=.3, ...)

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Matrix of dichotomous item responses. This matrix contain missing data (indicated
by \code{NA}) but missingness is assumed to be missing completely at random (MCAR).
}
  \item{weights}{
Optional vector of student weights.
}
  \item{Fval}{
Initial or fixed values of the loading matrix \eqn{\bold{F}}.
}
  \item{Fpatt}{
Pattern matrix of the loading matrix \eqn{\bold{F}}. If elements should be
estimated, then an entry of 1 must be included in the pattern matrix.
}
  \item{Pval}{
Initial or fixed values for the covariance matrix \eqn{\bold{P}}.
}
  \item{Ppatt}{
Pattern matrix for the covariance matrix \eqn{\bold{P}}.
}
  \item{Psival}{
Initial or fixed values for the matrix of residual correlations \eqn{\bold{\Psi}}.
}
  \item{Psipatt}{
Pattern matrix for the matrix of residual correlations \eqn{\bold{\Psi}}.
}
  \item{dimensions}{
Number of dimensions if an exploratory factor analysis should be estimated.
}
  \item{lower}{
Fixed vector of lower asymptotes \eqn{c_i}.
}
  \item{upper}{
Fixed vector of upper asymptotes \eqn{d_i}.
}
  \item{wgtm}{
Matrix with positive entries which indicates by a positive entry which
item pairs should be used for estimation.
}
\item{modesttype}{
Estimation type. \code{modesttype=1} refers to the NOHARM approximation,
while \code{modesttype=2} refers to the estimation based on tetrachoric
correlations.
}
\item{pos.loading}{An optional logical indicating whether all entries in the
loading matrix \eqn{\bold{F}} should be positive}
\item{pos.variance}{An optional logical indicating whether all variances (i.e.
    diagonal entries in \eqn{\bold{P}}) should be positive}
\item{pos.residcorr}{An optional logical indicating whether all
    entries in the matrix of residual correlations \eqn{\bold{\Psi}} should
    be positive}
  \item{maxiter}{
Maximum number of iterations
}
  \item{conv}{
Convergence criterion for parameters
}
  \item{increment.factor}{
Numeric value larger than 1 which controls the size of increments
in increasing iteration numbers. With a larger value, the increments
are heavier penalized.
}
\item{reliability}{Logical indicating whether reliability should be computed.}
\item{object}{Object of class \code{noharm.sirt}}
\item{file}{String indicating a file name for summary.}
\item{\dots}{Further arguments to be passed. }
}

%\item{x}{Object of class \code{noharm.sirt}}
%\item{what}{A character indicating whether the original estimates (\code{est})
%or standardized estimates (\code{stand}) should be displayed.}
%\item{efa.load.min}{Minimum absolute factor loading to be displayed in
%    plot for exploratory factor analysis.
%}

% For \code{plot.noharm.sirt},
% additional arguments for \code{\link[semPlot:semPaths]{semPlot::semPaths}} (from
%\pkg{semPlot} package) can be specified.

\details{
The NOHARM item response model follows the response equation
\deqn{P( X_{pi}=1 | \bold{\theta}_p )=c_i + ( d_i - c_i )
\Phi( f_{i0} + f_{i1} \theta_{p1} + ... +
f_{iD} \theta_{pD} ) }
for item responses \eqn{X_{pi}} of person \eqn{p} on
item \eqn{i}, \eqn{\bold{F}=(f_{id})} is a loading matrix and \eqn{\bold{P}}
the covariance matrix of \eqn{\bold{\theta}_p}. The lower
asymptotes \eqn{c_i} and upper asymptotes \eqn{d_i} must be
provided as fixed values.
The response equation can be equivalently written by introducing a latent
continuous item response \eqn{X_{pi}^\ast}
\deqn{ X_{pi}^\ast=f_{i0} + f_{i1} \theta_{p1} + ... +
f_{iD} \theta_{pD} + e_{pi} }
with a standard normally distributed residual \eqn{e_{pi}}. These residuals
have a correlation matrix \eqn{\bold{\Psi}} with ones in the diagonal.
In this \R implementation of the NOHARM model, correlations between residuals
are allowed.

See References for more details about estimation.
}


\value{
A list. The most important entries are
\item{tanaka}{Tanaka fit statistic}
\item{rmsr}{RMSR fit statistic}
\item{N.itempair}{Sample size per item pair}
\item{pm}{Product moment matrix}
\item{wgtm}{Matrix of weights for each item pair}
\item{sumwgtm}{Sum of lower triangle matrix \code{wgtm}}
\item{lower}{Lower asymptotes}
\item{upper}{Upper asymptotes}
\item{residuals}{Residual matrix from approximation of the \code{pm} matrix}
\item{final.constants}{Final constants}
\item{factor.cor}{Covariance matrix}
\item{thresholds}{Threshold parameters}
\item{uniquenesses}{Uniquenesses}
\item{loadings}{Matrix of standardized factor loadings (delta parametrization)}
\item{loadings.theta}{Matrix of factor loadings \eqn{\bold{F}} (theta parametrization)}
\item{residcorr}{Matrix of residual correlations}
\item{Nobs}{Number of observations}
\item{Nitems}{Number of items}
\item{Fpatt}{Pattern loading matrix for \eqn{\bold{F}}}
\item{Ppatt}{Pattern loading matrix for \eqn{\bold{P}}}
\item{Psipatt}{Pattern loading matrix for \eqn{\bold{\Psi}}}
\item{dat}{Used dataset}
\item{dimensions}{Number of dimensions}
\item{iter}{Number of iterations}
\item{Nestpars}{Number of estimated parameters}
\item{chisquare}{Statistic \eqn{\chi^2}}
\item{df}{Degrees of freedom}
\item{chisquare_df}{Ratio \eqn{\chi^2 / df}}
\item{rmsea}{RMSEA statistic}
\item{p.chisquare}{Significance for \eqn{\chi^2} statistic}
\item{omega.rel}{Reliability of the sum score according to Green and Yang (2009)}
}




\references{
Fraser, C., & McDonald, R. P. (1988). NOHARM: Least squares item factor analysis.
\emph{Multivariate Behavioral Research, 23}, 267-269.

Fraser, C., & McDonald, R. P. (2012). \emph{NOHARM 4 Manual}.
{http://noharm.niagararesearch.ca/nh4man/nhman.html}.

Green, S. B., & Yang, Y. (2009). Reliability of summed item
scores using structural equation modeling: An alternative to
coefficient alpha. \emph{Psychometrika, 74}, 155-167.

McDonald, R. P. (1982a). Linear versus models in item response theory.
\emph{Applied Psychological Measurement, 6}, 379-396.

McDonald, R. P. (1982b). \emph{Unidimensional and multidimensional models for
item response theory}. I.R.T., C.A.T. conference, Minneapolis, 1982, Proceedings.

McDonald, R. P. (1997). Normal-ogive multidimensional model.
In W. van der Linden & R. K. Hambleton (1997):
\emph{Handbook of modern item response theory} (pp. 257-269).
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
EAP person parameter estimates can be obtained by \code{\link{R2noharm.EAP}}.

Model fit can be assessed by \code{\link{modelfit.sirt}}.

See \code{\link{R2noharm}} for running the NOHARM software from within \R.

See Fraser and McDonald (2012) for an implementation of the NOHARM model which
is available as freeware ({http://noharm.niagararesearch.ca/};
the link seems to be broken in the meanwhile).
}


\examples{
#############################################################################
# EXAMPLE 1: Two-dimensional IRT model with 10 items
#############################################################################

#**** data simulation
set.seed(9776)
N <- 3400 # sample size
# define difficulties
f0 <- c( .5, .25, -.25, -.5, 0, -.5, -.25, .25, .5, 0 )
I <- length(f0)
# define loadings
f1 <- matrix( 0, I, 2 )
f1[ 1:5,1] <- c(.8,.7,.6,.5, .5)
f1[ 6:10,2] <- c(.8,.7,.6,.5, .5 )
# covariance matrix
Pval <- matrix( c(1,.5,.5,1), 2, 2 )
# simulate theta
library(mvtnorm)
theta <- mvtnorm::rmvnorm(N, mean=c(0,0), sigma=Pval )
# simulate item responses
dat <- matrix( NA, N, I )
for (ii in 1:I){ # ii <- 1
    dat[,ii] <- 1*( stats::pnorm(f0[ii]+theta[,1]*f1[ii,1]+theta[,2]*f1[ii,2])>
                     stats::runif(N) )
        }
colnames(dat) <- paste0("I", 1:I)

#**** Model 1: Two-dimensional CFA with estimated item loadings
# define pattern matrices
Pval <- .3+0*Pval
Ppatt <- 1*(Pval>0)
diag(Ppatt) <- 0
diag(Pval) <- 1
Fval <- .7 * ( f1>0)
Fpatt <- 1 * ( Fval > 0 )
# estimate model
mod1 <- sirt::noharm.sirt( dat=dat, Ppatt=Ppatt,Fpatt=Fpatt, Fval=Fval, Pval=Pval )
summary(mod1)
# EAP ability estimates
pmod1 <- sirt::R2noharm.EAP(mod1, theta.k=seq(-4,4,len=10) )
# model fit
summary( sirt::modelfit.sirt( mod1 ))

# estimate model based on tetrachoric correlations
mod1b <- sirt::noharm.sirt( dat=dat, Ppatt=Ppatt,Fpatt=Fpatt, Fval=Fval, Pval=Pval,
            modesttype=2)
summary(mod1b)

\dontrun{
#**** Model 2: Two-dimensional CFA with correlated residuals
# define pattern matrix for residual correlation
Psipatt <- 0*diag(I)
Psipatt[1,2] <- 1
Psival <- 0*Psipatt
# estimate model
mod2 <- sirt::noharm.sirt( dat=dat, Ppatt=Ppatt,Fpatt=Fpatt, Fval=Fval, Pval=Pval,
            Psival=Psival, Psipatt=Psipatt )
summary(mod2)

#**** Model 3: Two-dimensional Rasch model
# pattern matrices
Fval <- matrix(0,10,2)
Fval[1:5,1] <- Fval[6:10,2] <- 1
Fpatt <- 0*Fval
Ppatt <- Pval <- matrix(1,2,2)
Pval[1,2] <- Pval[2,1] <- 0
# estimate model
mod3 <- sirt::noharm.sirt( dat=dat, Ppatt=Ppatt,Fpatt=Fpatt, Fval=Fval, Pval=Pval )
summary(mod3)
# model fit
summary( sirt::modelfit.sirt( mod3 ))

#** compare fit with NOHARM
noharm.path <- "c:/NOHARM"
P.pattern <- Ppatt ; P.init <- Pval
F.pattern <- Fpatt ; F.init <- Fval
mod3b <- sirt::R2noharm( dat=dat, model.type="CFA",
             F.pattern=F.pattern, F.init=F.init, P.pattern=P.pattern,
             P.init=P.init, writename="example_sim_2dim_rasch",
             noharm.path=noharm.path, dec="," )
summary(mod3b)

#############################################################################
# EXAMPLE 2: data.read
#############################################################################

data(data.read)
dat <- data.read
I <- ncol(dat)

#**** Model 1: Unidimensional Rasch model
Fpatt <- matrix( 0, I, 1 )
Fval <- 1 + 0*Fpatt
Ppatt <- Pval <- matrix(1,1,1)
# estimate model
mod1 <- sirt::noharm.sirt( dat=dat, Ppatt=Ppatt,Fpatt=Fpatt, Fval=Fval, Pval=Pval )
summary(mod1)
plot(mod1)    # semPaths plot

#**** Model 2: Rasch model in which item pairs within a testlet are excluded
wgtm <- matrix( 1, I, I )
wgtm[1:4,1:4] <- wgtm[5:8,5:8] <- wgtm[ 9:12, 9:12] <- 0
# estimation
mod2 <- sirt::noharm.sirt(dat=dat, Ppatt=Ppatt,Fpatt=Fpatt, Fval=Fval, Pval=Pval, wgtm=wgtm)
summary(mod2)

#**** Model 3: Rasch model with correlated residuals
Psipatt <- Psival <- 0*diag(I)
Psipatt[1:4,1:4] <- Psipatt[5:8,5:8] <- Psipatt[ 9:12, 9:12] <- 1
diag(Psipatt) <- 0
Psival <- .6*(Psipatt>0)
# estimation
mod3 <- sirt::noharm.sirt( dat=dat, Ppatt=Ppatt,Fpatt=Fpatt, Fval=Fval, Pval=Pval,
            Psival=Psival, Psipatt=Psipatt )
summary(mod3)
# allow only positive residual correlations
mod3b <- sirt::noharm.sirt( dat=dat, Ppatt=Ppatt, Fpatt=Fpatt, Fval=Fval, Pval=Pval,
            Psival=Psival, Psipatt=Psipatt, pos.residcorr=TRUE)
summary(mod3b)

#**** Model 4: Rasch testlet model
Fval <- Fpatt <- matrix( 0, I, 4 )
Fval[,1] <- Fval[1:4,2] <- Fval[5:8,3] <- Fval[9:12,4 ] <- 1
Ppatt <- Pval <- diag(4)
colnames(Ppatt) <- c("g", "A", "B","C")
Pval <- .5*Pval
# estimation
mod4 <- sirt::noharm.sirt( dat=dat, Ppatt=Ppatt,Fpatt=Fpatt, Fval=Fval, Pval=Pval  )
summary(mod4)
# allow only positive variance entries
mod4b <- sirt::noharm.sirt( dat=dat, Ppatt=Ppatt,Fpatt=Fpatt, Fval=Fval, Pval=Pval,
               pos.variance=TRUE )
summary(mod4b)

#**** Model 5: Bifactor model
Fval <- matrix( 0, I, 4 )
Fval[,1] <- Fval[1:4,2] <- Fval[5:8,3] <- Fval[9:12,4 ] <- .6
Fpatt <- 1 * ( Fval > 0 )
Pval <- diag(4)
Ppatt <- 0*Pval
colnames(Ppatt) <- c("g", "A", "B","C")
# estimation
mod5 <- sirt::noharm.sirt( dat=dat, Ppatt=Ppatt,Fpatt=Fpatt, Fval=Fval, Pval=Pval  )
summary(mod5)
# allow only positive loadings
mod5b <- sirt::noharm.sirt( dat=dat, Ppatt=Ppatt,Fpatt=Fpatt, Fval=Fval, Pval=Pval,
              pos.loading=TRUE )
summary(mod5b)
summary( sirt::modelfit.sirt(mod5b))

#**** Model 6: 3-dimensional Rasch model
Fval <- matrix( 0, I, 3 )
Fval[1:4,1] <- Fval[5:8,2] <- Fval[9:12,3 ] <- 1
Fpatt <- 0*Fval
Pval <- .6*diag(3)
diag(Pval) <- 1
Ppatt <- 1+0*Pval
colnames(Ppatt) <- c("A", "B","C")
# estimation
mod6 <- sirt::noharm.sirt( dat=dat, Ppatt=Ppatt,Fpatt=Fpatt, Fval=Fval, Pval=Pval  )
summary(mod6)
summary( sirt::modelfit.sirt(mod6) )  # model fit

#**** Model 7: 3-dimensional 2PL model
Fval <- matrix( 0, I, 3 )
Fval[1:4,1] <- Fval[5:8,2] <- Fval[9:12,3 ] <- 1
Fpatt <- Fval
Pval <- .6*diag(3)
diag(Pval) <- 1
Ppatt <- 1+0*Pval
diag(Ppatt) <- 0
colnames(Ppatt) <- c("A", "B","C")
# estimation
mod7 <- sirt::noharm.sirt( dat=dat, Ppatt=Ppatt,Fpatt=Fpatt, Fval=Fval, Pval=Pval  )
summary(mod7)
summary( sirt::modelfit.sirt(mod7) )

#**** Model 8: Exploratory factor analysis with 3 dimensions
# estimation
mod8 <- sirt::noharm.sirt( dat=dat, dimensions=3  )
summary(mod8)
}
}

% plot(mod8)
% # plot only standardized loadings larger than .40
% plot(mod8, efa.load.min=.4 )

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{NOHARM}
\keyword{summary}
% \keyword{plot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
