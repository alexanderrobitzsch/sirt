%% File Name: rm.facets.Rd
%% File Version: 0.81

\name{rm.facets}
\alias{rm.facets}
\alias{summary.rm.facets}
\alias{logLik.rm.facets}
\alias{anova.rm.facets}
\alias{IRT.irfprob.rm.facets}
\alias{IRT.likelihood.rm.facets}
\alias{IRT.posterior.rm.facets}
\alias{IRT.modelfit.rm.facets}
\alias{IRT.factor.scores.rm.facets}
\alias{summary.IRT.modelfit.rm.facets}
\alias{rm_proc_data}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rater Facets Models with Item/Rater Intercepts and Slopes
}
\description{
This function estimates the unidimensional rater facets model (Lincare, 1994)
and an extension to slopes (see Details; Robitzsch & Steinfeld, 2018). The estimation
is conducted by an EM algorithm employing marginal
maximum likelihood.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\usage{
rm.facets(dat, pid=NULL, rater=NULL, Qmatrix=NULL, theta.k=seq(-9, 9, len=30),
    est.b.rater=TRUE, est.a.item=FALSE, est.a.rater=FALSE, rater_item_int=FALSE,
    est.mean=FALSE, tau.item.fixed=NULL, a.item.fixed=NULL, b.rater.fixed=NULL,
    a.rater.fixed=NULL, b.rater.center=2, a.rater.center=2, a.item.center=2, a_lower=.05,
    a_upper=10, reference_rater=NULL, max.b.increment=1, numdiff.parm=0.00001,
    maxdevchange=0.1, globconv=0.001, maxiter=1000, msteps=4, mstepconv=0.001,
    PEM=FALSE, PEM_itermax=maxiter)

\method{summary}{rm.facets}(object, file=NULL, ...)

\method{anova}{rm.facets}(object,...)

\method{logLik}{rm.facets}(object,...)

\method{IRT.irfprob}{rm.facets}(object,...)

\method{IRT.factor.scores}{rm.facets}(object, type="EAP", ...)

\method{IRT.likelihood}{rm.facets}(object,...)

\method{IRT.posterior}{rm.facets}(object,...)

\method{IRT.modelfit}{rm.facets}(object,...)

\method{summary}{IRT.modelfit.rm.facets}(object, ...)

## function for processing data
rm_proc_data( dat, pid, rater, rater_item_int=FALSE, reference_rater=NULL )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Original data frame. Ratings on variables must be in rows,
i.e. every row corresponds to a person-rater combination.
}
  \item{pid}{
Person identifier.
}
  \item{rater}{
Rater identifier
}
  \item{Qmatrix}{
An optional Q-matrix. If this matrix is not provided,
then by default the ordinary scoring of categories
(from 0 to the maximum score of \eqn{K}) is used.
}
  \item{theta.k}{
A grid of theta values for the ability distribution.
}
  \item{est.b.rater}{
Should the rater severities \eqn{b_r} be estimated?
}
  \item{est.a.item}{
Should the item slopes \eqn{a_i} be estimated?
}
  \item{est.a.rater}{
Should the rater slopes \eqn{a_r} be estimated?
}
\item{rater_item_int}{Logical indicating whether rater-item-interactions
should be modeled.}
\item{est.mean}{Optional logical indicating whether the mean of the
trait distribution should be estimated.}
\item{tau.item.fixed}{Matrix with fixed \eqn{\tau} parameters. Non-fixed
parameters must be declared by \code{NA} values.}
\item{a.item.fixed}{Vector with fixed item discriminations}
\item{b.rater.fixed}{Vector with fixed rater intercept parameters}
\item{a.rater.fixed}{Vector with fixed rater discrimination parameters}
\item{b.rater.center}{Centering method for rater intercept parameters. The
    value \code{0} corresponds to no centering, the values \code{1} and
\code{2} to different methods to ensure that they sum to zero.}
\item{a.rater.center}{Centering method for rater discrimination parameters. The
    value \code{0} corresponds to no centering, the values \code{1} and
\code{2} to different methods to ensure that their product equals one.}
\item{a.item.center}{Centering method for item discrimination parameters. The
    value \code{0} corresponds to no centering, the values \code{1} and
\code{2} to different methods to ensure that their product equals one.}
\item{a_lower}{Lower bound for \eqn{a} parameters}
\item{a_upper}{Upper bound for \eqn{a} parameters}
\item{reference_rater}{Identifier for rater as a reference rater for which
a fixed rater mean of 0 and a fixed rater slope of 1 is assumed.}
\item{max.b.increment}{
Maximum increment of item parameters during estimation
}
  \item{numdiff.parm}{
Numerical differentiation step width
}
  \item{maxdevchange}{
Maximum relative deviance change as a convergence criterion
}
  \item{globconv}{
Maximum parameter change
}
  \item{maxiter}{
Maximum number of iterations
}
  \item{msteps}{
Maximum number of iterations during an M step
}
  \item{mstepconv}{
Convergence criterion in an M step
}
\item{PEM}{Logical indicating whether the P-EM acceleration should be
applied (Berlinet & Roland, 2012).}
\item{PEM_itermax}{Number of iterations in which the P-EM method should be
applied.}
  \item{object}{
Object of class \code{rm.facets}
}
\item{file}{Optional file name in which summary should be written.}
\item{type}{Factor score estimation method. Factor score types
    \code{"EAP"}, \code{"MLE"} and \code{"WLE"} are supported.}
  \item{\dots}{
Further arguments to be passed
}
}


% \item{skillspace}{Type of skill space which is used to smooth the distribution
%    on nodes \code{theta.k}. \code{skillspace="normal"} corresponds to
%    a normal distribution, \code{skillspace="loglinear"} to log-linear
%    modelling of the skill space distribution for which a design matrix \eqn{Z}
%    has to be specified, and \code{skillspace="no"} corresponds to an
% unsmoothed distribution.}

%\item{Z}{Design matrix for trait distribution if \code{skillspace="loglinear"}.
%By using appropriate design matrices, non-normal trait distributions
%can be estimated.}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETAILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\details{
This function models ratings \eqn{X_{pri}}
for person \eqn{p}, rater \eqn{r} and item \eqn{i}
and category \eqn{k} (see also Robitzsch & Steinfeld, 2018; Uto & Ueno, 2010; Wu, 2017)
\deqn{P( X_{pri}=k | \theta_p ) \propto
    exp( a_i a_r q_{ik} \theta_p - q_{ik} b_r -   \tau_{ik} ) \quad,
    \quad \theta_p \sim N( 0, \sigma^2 )}{%
    P( X_{pri}=k | \theta_p ) \propto
    exp( a_i a_r q_{ik} \theta_p - q_{ik} b_r - \tau_{ik} ),
            \theta_p ~ N( 0, \sigma^2 )}
By default, the scores in the \eqn{Q} matrix are
\eqn{q_{ik}=k}. Item slopes \eqn{a_i} and rater slopes
\eqn{a_r} are standardized such that their product equals
one, i.e. \eqn{ \prod_i a_i=\prod_r a_r=1}.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\value{
A list with following entries:
\item{deviance}{Deviance}
\item{ic}{Information criteria and number of parameters}
\item{item}{Data frame with item parameters}
\item{rater}{Data frame with rater parameters}
\item{person}{Data frame with person parameters: EAP and corresponding
    standard errors}
\item{EAP.rel}{EAP reliability}
\item{mu}{Mean of the trait distribution}
\item{sigma}{Standard deviation of the trait distribution}
\item{theta.k}{Grid of theta values}
\item{pi.k}{Fitted distribution at \code{theta.k} values}
\item{tau.item}{Item parameters \eqn{\tau_{ik}}}
\item{se.tau.item}{Standard error of item parameters \eqn{\tau_{ik}}}
\item{a.item}{Item slopes \eqn{a_i}}
\item{se.a.item}{Standard error of item slopes \eqn{a_i}}
\item{delta.item}{Delta item parameter. See
    \code{\link{pcm.conversion}}.
       }
\item{b.rater}{Rater severity parameter \eqn{b_r}}
\item{se.b.rater}{Standard error of rater severity parameter \eqn{b_r}}
\item{a.rater}{Rater slope parameter \eqn{a_r}}
\item{se.a.rater}{Standard error of rater slope parameter \eqn{a_r}}
\item{f.yi.qk}{Individual likelihood}
\item{f.qk.yi}{Individual posterior distribution}
\item{probs}{Item probabilities at grid \code{theta.k}}
\item{n.ik}{Expected counts}
\item{maxK}{Maximum number of categories}
\item{procdata}{Processed data}
\item{iter}{Number of iterations}
\item{ipars.dat2}{Item parameters for expanded dataset \code{dat2}}
\item{\dots}{Further values}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\references{
Berlinet, A. F., & Roland, C. (2012).
Acceleration of the EM algorithm: P-EM versus epsilon algorithm.
\emph{Computational Statistics & Data Analysis, 56}(12), 4122-4137.

Linacre, J. M. (1994). \emph{Many-Facet Rasch Measurement}.
Chicago: MESA Press.

Robitzsch, A., & Steinfeld, J. (2018). Item response models for human ratings: Overview,
estimation methods, and implementation in R.
\emph{Psychological Test and Assessment Modeling, 60}(1), 101-139.

Uto, M., & Ueno, M. (2016). Item response theory for peer assessment.
\emph{IEEE Transactions on Learning Technologies, 9}(2), 157-170.

Wu, M. (2017). Some IRT-based analyses for interpreting rater effects.
\emph{Psychological Test and Assessment Modeling, 59}(4), 453-470.
}


\author{
Alexander Robitzsch
}

\note{
If the trait standard deviation \code{sigma} strongly
differs from 1, then a user should investigate the sensitivity
of results using different theta integration points \code{theta.k}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also the \pkg{TAM} package for the estimation
of more complicated facet models.

See \code{\link{rm.sdt}} for estimating a hierarchical rater model.
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\examples{
#############################################################################
# EXAMPLE 1: Partial Credit Model and Generalized partial credit model
#                   5 items and 1 rater
#############################################################################
data(data.ratings1)
dat <- data.ratings1

# select rater db01
dat <- dat[ paste(dat$rater)=="db01", ]

#****  Model 1: Partial Credit Model
mod1 <- sirt::rm.facets( dat[, paste0( "k",1:5) ], pid=dat$idstud )

#****  Model 2: Generalized Partial Credit Model
mod2 <- sirt::rm.facets( dat[, paste0( "k",1:5) ],  pid=dat$idstud, est.a.item=TRUE)

summary(mod1)
summary(mod2)

\dontrun{
#############################################################################
# EXAMPLE 2: Facets Model: 5 items, 7 raters
#############################################################################

data(data.ratings1)
dat <- data.ratings1

#****  Model 1: Partial Credit Model: no rater effects
mod1 <- sirt::rm.facets( dat[, paste0( "k",1:5) ], rater=dat$rater,
             pid=dat$idstud, est.b.rater=FALSE )

#****  Model 2: Partial Credit Model: intercept rater effects
mod2 <- sirt::rm.facets( dat[, paste0( "k",1:5) ], rater=dat$rater, pid=dat$idstud)

# extract individual likelihood
lmod1 <- IRT.likelihood(mod1)
str(lmod1)
# likelihood value
logLik(mod1)
# extract item response functions
pmod1 <- IRT.irfprob(mod1)
str(pmod1)
# model comparison
anova(mod1,mod2)
# absolute and relative model fit
smod1 <- IRT.modelfit(mod1)
summary(smod1)
smod2 <- IRT.modelfit(mod2)
summary(smod2)
IRT.compareModels( smod1, smod2 )
# extract factor scores (EAP is the default)
IRT.factor.scores(mod2)
# extract WLEs
IRT.factor.scores(mod2, type="WLE")

#****  Model 2a: compare results with TAM package
#   Results should be similar to Model 2
library(TAM)
mod2a <- TAM::tam.mml.mfr( resp=dat[, paste0( "k",1:5) ],
             facets=dat[, "rater", drop=FALSE],
             pid=dat$pid, formulaA=~ item*step + rater )

#****  Model 2b: Partial Credit Model: some fixed parameters
# fix rater parameters for raters 1, 4 and 5
b.rater.fixed <- rep(NA,7)
b.rater.fixed[ c(1,4,5) ] <- c(1,-.8,0)  # fixed parameters
# fix item parameters of first and second item
tau.item.fixed <- round( mod2$tau.item, 1 )    # use parameters from mod2
tau.item.fixed[ 3:5, ] <- NA    # free item parameters of items 3, 4 and 5
mod2b <- sirt::rm.facets( dat[, paste0( "k",1:5) ], rater=dat$rater,
             b.rater.fixed=b.rater.fixed, tau.item.fixed=tau.item.fixed,
             est.mean=TRUE, pid=dat$idstud)
summary(mod2b)

#****  Model 3: estimated rater slopes
mod3 <- sirt::rm.facets( dat[, paste0( "k",1:5) ], rater=dat$rater,
            est.a.rater=TRUE)

#****  Model 4: estimated item slopes
mod4 <- sirt::rm.facets( dat[, paste0( "k",1:5) ], rater=dat$rater,
             pid=dat$idstud, est.a.item=TRUE)

#****  Model 5: estimated rater and item slopes
mod5 <- sirt::rm.facets( dat[, paste0( "k",1:5) ], rater=dat$rater,
             pid=dat$idstud, est.a.rater=TRUE, est.a.item=TRUE)
summary(mod1)
summary(mod2)
summary(mod2a)
summary(mod3)
summary(mod4)
summary(mod5)

#****  Model 5a: Some fixed parameters in Model 5
# fix rater b parameters for raters 1, 4 and 5
b.rater.fixed <- rep(NA,7)
b.rater.fixed[ c(1,4,5) ] <- c(1,-.8,0)
# fix rater a parameters for first four raters
a.rater.fixed <- rep(NA,7)
a.rater.fixed[ c(1,2,3,4) ] <- c(1.1,0.9,.85,1)
# fix item b parameters of first item
tau.item.fixed <- matrix( NA, nrow=5, ncol=3 )
tau.item.fixed[ 1, ] <- c(-2,-1.5, 1 )
# fix item a parameters
a.item.fixed <- rep(NA,5)
a.item.fixed[ 1:4 ] <- 1
# estimate model
mod5a <- sirt::rm.facets( dat[, paste0( "k",1:5) ], rater=dat$rater,
             pid=dat$idstud, est.a.rater=TRUE, est.a.item=TRUE,
             tau.item.fixed=tau.item.fixed, b.rater.fixed=b.rater.fixed,
             a.rater.fixed=a.rater.fixed, a.item.fixed=a.item.fixed,
             est.mean=TRUE)
summary(mod5a)

#****  Model 6: Estimate rater model with reference rater 'db03'
mod6 <- sirt::rm.facets( dat[, paste0( "k",1:5) ], rater=dat$rater, est.a.item=TRUE,
             est.a.rater=TRUE, pid=dat$idstud, reference_rater="db03" )
summary(mod6)

#**** Model 7: Modelling rater-item-interactions
mod7 <- sirt::rm.facets( dat[, paste0( "k",1:5) ], rater=dat$rater, est.a.item=FALSE,
             est.a.rater=TRUE, pid=dat$idstud, reference_rater="db03",
             rater_item_int=TRUE)
summary(mod7)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Rater model}
\keyword{Facets model}
\keyword{summary}
