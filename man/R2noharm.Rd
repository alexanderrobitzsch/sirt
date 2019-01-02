%% File Name: R2noharm.Rd
%% File Version: 2.33

\name{R2noharm}
\alias{R2noharm}
\alias{summary.R2noharm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of a NOHARM Analysis from within \R
}
\description{
This function enables the estimation of a NOHARM analysis
(Fraser & McDonald, 1988; McDonald, 1982a, 1982b, 1997) from within \R.
NOHARM estimates a compensatory multidimensional factor analysis for dichotomous
response data. Arguments of this function strictly follow
the rules of the NOHARM manual (see Fraser & McDonald, 2012).
}


\usage{
R2noharm(dat=NULL,pm=NULL, n=NULL, model.type, weights=NULL, dimensions=NULL,
      guesses=NULL, noharm.path, F.pattern=NULL, F.init=NULL,
      P.pattern=NULL, P.init=NULL, digits.pm=4, writename=NULL,
      display.fit=5,  dec=".", display=TRUE)

\method{summary}{R2noharm}(object, logfile=NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\arguments{
  \item{dat}{
An \eqn{N \times I} data frame of item responses for \eqn{N} subjects
and \eqn{I} items
}
\item{pm}{A matrix or a vector containing product-moment correlations}
\item{n}{Sample size. This value must only be included if \code{pm} is
provided.}
  \item{model.type}{
Can be \code{"EFA"} (exploratory factor analysis) or \code{"CFA"}
(confirmatory factor analysis).
}
\item{weights}{Optional vector of student weights}
  \item{dimensions}{
Number of dimensions in exploratory factor analysis
}
  \item{guesses}{
An optional vector of fixed guessing parameters of length \eqn{I}.
In case of the default \code{NULL}, all guessing parameters are set
to zero.
}
  \item{noharm.path}{
Local path where the NOHARM 4 command line 64-bit version is located.
}
  \item{F.pattern}{
Pattern matrix for \eqn{F} (\eqn{I \times D})
}
  \item{F.init}{
Initial matrix for \eqn{F} (\eqn{I \times D})
}
  \item{P.pattern}{
Pattern matrix for \eqn{P} (\eqn{D \times D})
}
  \item{P.init}{
Initial matrix for \eqn{P} (\eqn{D \times D})
}
  \item{digits.pm}{
Number of digits after decimal separator which are used for
estimation
}
  \item{writename}{
Name for NOHARM input and output files
  }
  \item{display.fit}{
How many digits (after decimal separator) should be used
for printing results on the \R console?
}
  \item{dec}{
Decimal separator (\code{"."} or \code{","})
}
  \item{display}{
Display output?
}
  \item{object}{
Object of class \code{R2noharm}
}
\item{logfile}{File name if the summary should be sunk into a file}
\item{\dots}{Further arguments to be passed}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETAILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\details{
NOHARM estimates a multidimensional compensatory
item response model with the probit link function \eqn{\Phi}.
For item responses \eqn{X_{pi}} of person \eqn{p} on
item \eqn{i} the model equation is defined as
\deqn{P( X_{pi}=1 | \bold{\theta}_p )=c_i + ( 1 - c_i )
\Phi( f_{i0} + f_{i1} \theta_{p1} + ... +
f_{iD} \theta_{pD} ) }
where \eqn{F=(f_{id})} is a loading matrix and \eqn{P}
the covariance matrix of \eqn{\bold{\theta}_p}. The guessing
parameters \eqn{c_i} must be provided as fixed values.

For the definition of \eqn{F} and \eqn{P} matrices, please
consult the NOHARM manual.

This function needs the 64-bit command line version which can be downloaded
at \emph{http://noharm.niagararesearch.ca/nh4cldl.html}.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\value{
A list with following entries
\item{tanaka}{Tanaka index}
\item{rmsr}{RMSR statistic}
\item{N.itempair}{Sample sizes of pairwise item observations}
\item{pm}{Product moment matrix}
\item{weights}{Used student weights}
\item{guesses}{Fixed guessing parameters}
\item{residuals}{Residual covariance matrix}
\item{final.constants}{Vector of final constants}
\item{thresholds}{Threshold parameters}
\item{uniquenesses}{Item uniquenesses}
\item{loadings.theta}{Matrix of loadings in theta parametrization
    (common factor parametrization)}
\item{factor.cor}{Covariance matrix of factors}
\item{difficulties}{Item difficulties (for unidimensional models)}
\item{discriminations}{Item discriminations (for unidimensional models)}
\item{loadings}{Loading matrix (latent trait parametrization)}
\item{model.type}{Used model type}
\item{Nobs}{Number of observations}
\item{Nitems}{Number of items}
\item{modtype}{Model type according to the NOHARM specification (see NOHARM manual)}
\item{F.init}{Initial loading matrix for \eqn{F}}
\item{F.pattern}{Pattern loading matrix for \eqn{F}}
\item{P.init}{Initial covariance matrix for \eqn{P}}
\item{P.pattern}{Pattern covariance matrix for \eqn{P}}
\item{dat}{Original data frame}
\item{systime}{System time}
\item{noharm.path}{Used NOHARM directory}
\item{digits.pm}{Number of digits in product moment matrix}
\item{dec}{Used decimal symbol}
\item{display.fit}{Number of digits for fit display}
\item{dimensions}{Number of dimensions}
\item{chisquare}{Statistic \eqn{\chi^2}}
\item{Nestpars}{Number of estimated parameters}
\item{df}{Degrees of freedom}
\item{chisquare_df}{Ratio \eqn{\chi^2 / df}}
\item{rmsea}{RMSEA statistic}
\item{p.chisquare}{Significance for \eqn{\chi^2} statistic}
}
\references{
Fraser, C., & McDonald, R. P. (1988). NOHARM: Least squares item factor analysis.
\emph{Multivariate Behavioral Research, 23}, 267-269.

Fraser, C., & McDonald, R. P. (2012). \emph{NOHARM 4 Manual}.
http://noharm.niagararesearch.ca/nh4man/nhman.html.

McDonald, R. P. (1982a). Linear versus models in item response theory.
\emph{Applied Psychological Measurement, 6}, 379-396.

McDonald, R. P. (1982b). \emph{Unidimensional and multidimensional models for
item response theory}. I.R.T., C.A.T. conference, Minneapolis, 1982, Proceedings.

McDonald, R. P. (1997). Normal-ogive multidimensional model.
In W. van der Linden & R. K. Hambleton (1997):
\emph{Handbook of modern item response theory} (pp. 257-269).
New York: Springer.
}



%\author{
%Alexander Robitzsch
%}


\note{
Possible errors often occur due to wrong \code{dec}
specification.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For estimating standard errors see \code{\link{R2noharm.jackknife}}.

For EAP person parameter estimates see \code{\link{R2noharm.EAP}}.

For an \R implementation of the NOHARM model see \code{\link{noharm.sirt}}.
}
%% ~~objects to See Also as \code{\link{help}}, ~~~


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Data data.noharm18 with 18 items
#############################################################################

# load data
data( data.noharm18 )
dat <- data.noharm18
I <- ncol(dat) # number of items

# locate noharm.path
noharm.path <- "c:/NOHARM"

#****************************************
# Model 1: 1-dimensional Rasch model (1-PL model)

# estimate one factor variance
P.pattern <- matrix( 1, ncol=1, nrow=1 )
P.init <- P.pattern
# fix all entries in the loading matrix to 1
F.pattern <- matrix( 0, nrow=I, ncol=1 )
F.init <- 1 + 0*F.pattern       #
# estimate model
mod <- sirt::R2noharm( dat=dat, model.type="CFA",
           F.pattern=F.pattern, F.init=F.init, P.pattern=P.pattern,
           P.init=P.init, writename="ex1__1dim_1pl",
       noharm.path=noharm.path, dec="," )
# summary
summary(mod, logfile="ex1__1dim_1pl__SUMMARY")
# jackknife
jmod <- sirt::R2noharm.jackknife( mod, jackunits=20 )
summary(jmod, logfile="ex1__1dim_1pl__JACKKNIFE")
# compute factor scores (EAPs)
emod <- sirt::R2noharm.EAP(mod)

#*****-----
# Model 1b: Include student weights in estimation
N <- nrow(dat)
weights <- stats::runif( N, 1, 5 )
mod1b <- sirt::R2noharm( dat=dat, model.type="CFA",  weights=weights,
            F.pattern=F.pattern, F.init=F.init, P.pattern=P.pattern,
            P.init=P.init, writename="ex1__1dim_1pl_w",
            noharm.path=noharm.path, dec="," )
summary(mod1b)

#****************************************
# Model 2: 1-dimensional 2PL Model

# set trait variance equal to 1
P.pattern <- matrix( 0, ncol=1, nrow=1 )
P.init <- 1+0*P.pattern
# loading matrix
F.pattern <- matrix( 1, nrow=I, ncol=1 )
F.init <- 1 + 0*F.pattern

mod <- sirt::R2noharm( dat=dat, model.type="CFA",
            F.pattern=F.pattern, F.init=F.init, P.pattern=P.pattern,
            P.init=P.init, writename="ex2__1dim_2pl",
            noharm.path=noharm.path, dec="," )

summary(mod)
jmod <- sirt::R2noharm.jackknife( mod, jackunits=20 )
summary(jmod)

#****************************************
# Model 3: 1-dimensional 3PL Model with fixed guessing parameters

# set trait variance equal to 1
P.pattern <- matrix( 0, ncol=1, nrow=1 )
P.init <- 1+0*P.pattern
# loading matrix
F.pattern <- matrix( 1, nrow=I, ncol=1 )
F.init <- 1 + 0*F.pattern       #
# fix guessing parameters equal to .2 (for all items)
guesses <- rep( .1, I )

mod <- sirt::R2noharm( dat=dat, model.type="CFA",
          F.pattern=F.pattern, F.init=F.init, P.pattern=P.pattern,
          P.init=P.init, guesses=guesses,
          writename="ex3__1dim_3pl", noharm.path=noharm.path, dec=","  )
summary(mod)
jmod <- sirt::R2noharm.jackknife( mod, jackunits=20 )
summary(jmod)

#****************************************
# Model 4: 3-dimensional Rasch model

# estimate one factor variance
P.pattern <- matrix( 1, ncol=3, nrow=3 )
P.init <- .8*P.pattern
diag(P.init) <- 1
# fix all entries in the loading matrix to 1
F.init <- F.pattern <- matrix( 0, nrow=I, ncol=3 )
F.init[1:6,1] <- 1
F.init[7:12,2] <- 1
F.init[13:18,3] <- 1

mod <- sirt::R2noharm( dat=dat, model.type="CFA",
          F.pattern=F.pattern, F.init=F.init, P.pattern=P.pattern,
          P.init=P.init, writename="ex4__3dim_1pl",
          noharm.path=noharm.path, dec="," )
# write output from R console in a file
summary(mod, logfile="ex4__3dim_1pl__SUMMARY.Rout")

jmod <- sirt::R2noharm.jackknife( mod, jackunits=20 )
summary(jmod)

# extract factor scores
emod <- sirt::R2noharm.EAP(mod)

#****************************************
# Model 5: 3-dimensional 2PL model

# estimate one factor variance
P.pattern <- matrix( 1, ncol=3, nrow=3 )
P.init <- .8*P.pattern
diag(P.init) <- 0
# fix all entries in the loading matrix to 1
F.pattern <- matrix( 0, nrow=I, ncol=3 )
F.pattern[1:6,1] <- 1
F.pattern[7:12,2] <- 1
F.pattern[13:18,3] <- 1
F.init <- F.pattern

mod <- sirt::R2noharm( dat=dat, model.type="CFA",
          F.pattern=F.pattern, F.init=F.init, P.pattern=P.pattern,
          P.init=P.init, writename="ex5__3dim_2pl",
          noharm.path=noharm.path, dec="," )
summary(mod)
# use 50 jackknife units with 4 persons within a unit
jmod <- sirt::R2noharm.jackknife( mod, jackunits=seq( 1:50, each=4 ) )
summary(jmod)

#****************************************
# Model 6: Exploratory Factor Analysis with 3 factors

mod <- sirt::R2noharm( dat=dat, model.type="EFA",  dimensions=3,
           writename="ex6__3dim_efa", noharm.path=noharm.path,dec=",")
summary(mod)

jmod <- sirt::R2noharm.jackknife( mod, jackunits=20 )

#############################################################################
# EXAMPLE 2: NOHARM manual Example A
#############################################################################

# See NOHARM manual: http://noharm.niagararesearch.ca/nh4man/nhman.html
# The following text and data is copied from this manual.
#
# In the first example, we demonstrate how to prepare the input for a 2-dimensional
# model using exploratory analysis. Data from a 9 item test were collected from
# 200 students and the 9x9 product-moment matrix of the responses was computed.
#
# Our hypothesis is for a 2-dimensional model with no guessing,
# i.e., all guesses are equal to zero. However, because we are unsure of any
# particular pattern for matrix F, we wish to prescribe an exploratory analysis, i.e.,
# set EX=1. Also, we will content ourselves with letting the program supply all
# initial values.
#
# We would like both the sample product-moment matrix and the residual matrix to
# be included in the output.

# scan product-moment matrix copied from the NOHARM manual
pm <- scan()
     0.8967
     0.2278 0.2356
     0.6857 0.2061 0.7459
     0.8146 0.2310 0.6873 0.8905
     0.4505 0.1147 0.3729 0.4443 0.5000
     0.7860 0.2080 0.6542 0.7791 0.4624 0.8723
     0.2614 0.0612 0.2140 0.2554 0.1914 0.2800 0.2907
     0.7549 0.1878 0.6236 0.7465 0.4505 0.7590 0.2756 0.8442
     0.6191 0.1588 0.5131 0.6116 0.3845 0.6302 0.2454 0.6129 0.6879

ex2 <- sirt::R2noharm( pm=pm, n=200, model.type="EFA", dimensions=2,
         noharm.path=noharm.path, writename="ex2_noharmExA", dec=",")
summary(ex2)

#############################################################################
# EXAMPLE 3: NOHARM manual Example B
#############################################################################

# See NOHARM manual: http://noharm.niagararesearch.ca/nh4man/nhman.html
# The following text and data is copied from this manual.

# Suppose we have the product-moment matrix of data from 125 students on 9 items.
# Our hypothesis is for 2 dimensions with simple structure. In this case,
# items 1 to 5 have coefficients of theta which are to be estimated for one
# latent trait but are to be fixed at zero for the other one.
# For the latent trait for which items 1 to 5 have zero coefficients,
# items 6 to 9 have coefficients which are to be estimated. For the other
# latent trait, items 6 to 9 will have zero coefficients.
# We also wish to estimate the correlation between the latent traits,
# so we prescribe P as a 2x2 correlation matrix.
#
# Our hypothesis prescribes that there was no guessing involved, i.e.,
# all guesses are equal to zero. For demonstration purposes,
# let us not have the program print out the sample product-moment matrix.
# Also let us not supply any starting values but, rather, use the defaults
# supplied by the program.

pm <- scan()
    0.930
    0.762 0.797
    0.541 0.496 0.560
    0.352 0.321 0.261 0.366
    0.205 0.181 0.149 0.110 0.214
    0.858 0.747 0.521 0.336 0.203 0.918
    0.773 0.667 0.465 0.308 0.184 0.775 0.820
    0.547 0.474 0.347 0.233 0.132 0.563 0.524 0.579
    0.329 0.290 0.190 0.140 0.087 0.333 0.308 0.252 0.348

I <- 9    # number of items
# define loading matrix
F.pattern <- matrix(0,I,2)
F.pattern[1:5,1] <- 1
F.pattern[6:9,2] <- 1
F.init <- F.pattern
# define covariance matrix
P.pattern <- matrix(1,2,2)
diag(P.pattern) <- 0
P.init <- 1+P.pattern

ex3 <- sirt::R2noharm( pm=pm, n=125,, model.type="CFA",
           F.pattern=F.pattern, F.init=F.init, P.pattern=P.pattern,
           P.init=P.init, writename="ex3_noharmExB",
           noharm.path=noharm.path, dec="," )
summary(ex3)

#############################################################################
# EXAMPLE 4: NOHARM manual Example C
#############################################################################

data(data.noharmExC)
# See NOHARM manual: http://noharm.niagararesearch.ca/nh4man/nhman.html
# The following text and data is copied from this manual.

# In this example, suppose that from 300 respondents we have item
# responses scored dichotomously, 1 or 0, for 8 items.
#
# Our hypothesis is for a unidimensional model where all eight items
# have coefficients of theta which are to be estimated.
# Suppose that since the items were multiple choice with 5 options each,
# we set the fixed guesses all to 0.2 (not necessarily good reasoning!)
#
# Let's supply initial values for the coefficients of theta (F matrix)
# as .75 for items 1 to 4 and .6 for items 5 to 8.

I <- 8
guesses <- rep(.2,I)
F.pattern <- matrix(1,I,1)
F.init <- F.pattern
F.init[1:4,1] <- .75
F.init[5:8,1] <- .6
P.pattern <- matrix(0,1,1)
P.init <- 1 + 0 * P.pattern

ex4 <- sirt::R2noharm( dat=data.noharmExC,, model.type="CFA",
           guesses=guesses, F.pattern=F.pattern, F.init=F.init,
           P.pattern=P.pattern, P.init=P.init, writename="ex3_noharmExC",
           noharm.path=noharm.path, dec="," )
summary(ex4)

# modify F pattern matrix
# f11=f51 (since both have equal pattern values of 2),
# f21=f61 (since both have equal pattern values of 3),
# f31=f71 (since both have equal pattern values of 4),
# f41=f81 (since both have equal pattern values of 5).
F.pattern[ c(1,5) ] <- 2
F.pattern[ c(2,6) ] <- 3
F.pattern[ c(3,7) ] <- 4
F.pattern[ c(4,8) ] <- 5
F.init <- .5+0*F.init

ex4a <- sirt::R2noharm( dat=data.noharmExC,, model.type="CFA",
           guesses=guesses, F.pattern=F.pattern, F.init=F.init,
           P.pattern=P.pattern, P.init=P.init, writename="ex3_noharmExC1",
           noharm.path=noharm.path, dec="," )
summary(ex4a)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{NOHARM}
%% \keyword{summary}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
