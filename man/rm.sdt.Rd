%% File Name: rm.sdt.Rd
%% File Version: 1.48

\name{rm.sdt}
\alias{rm.sdt}
\alias{summary.rm.sdt}
\alias{logLik.rm.sdt}
\alias{anova.rm.sdt}
\alias{IRT.irfprob.rm.sdt}
\alias{IRT.likelihood.rm.sdt}
\alias{IRT.posterior.rm.sdt}
\alias{IRT.modelfit.rm.sdt}
\alias{summary.IRT.modelfit.rm.sdt}
\alias{IRT.factor.scores.rm.sdt}
\alias{plot.rm.sdt}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hierarchical Rater Model Based on Signal Detection Theory (HRM-SDT)
}

\description{
This function estimates a version of the hierarchical rater
model (HRM) based on signal detection theory (HRM-SDT; DeCarlo, 2005;
DeCarlo, Kim & Johnson, 2011; Robitzsch & Steinfeld, 2018).
The model is estimated by means of an EM algorithm adapted from
multilevel latent class analysis (Vermunt, 2008).
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\usage{
rm.sdt(dat, pid, rater, Qmatrix=NULL, theta.k=seq(-9, 9, len=30),
    est.a.item=FALSE, est.c.rater="n", est.d.rater="n", est.mean=FALSE, est.sigma=TRUE,
    skillspace="normal", tau.item.fixed=NULL, a.item.fixed=NULL,
    d.min=0.5, d.max=100, d.start=3, c.start=NULL, tau.start=NULL, sd.start=1,
    d.prior=c(3,100), c.prior=c(3,100), tau.prior=c(0,1000), a.prior=c(1,100),
    link_item="GPCM", max.increment=1, numdiff.parm=0.00001, maxdevchange=0.1,
    globconv=.001, maxiter=1000, msteps=4, mstepconv=0.001, optimizer="nlminb" )

\method{summary}{rm.sdt}(object, file=NULL, ...)

\method{plot}{rm.sdt}(x, ask=TRUE, ...)

\method{anova}{rm.sdt}(object,...)

\method{logLik}{rm.sdt}(object,...)

\method{IRT.factor.scores}{rm.sdt}(object, type="EAP", ...)

\method{IRT.irfprob}{rm.sdt}(object,...)

\method{IRT.likelihood}{rm.sdt}(object,...)

\method{IRT.posterior}{rm.sdt}(object,...)

\method{IRT.modelfit}{rm.sdt}(object,...)

\method{summary}{IRT.modelfit.rm.sdt}(object,...)
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
Rater identifier.
}
  \item{Qmatrix}{
An optional Q-matrix. If this matrix is not provided,
then by default the ordinary scoring of categories
(from 0 to the maximum score of \eqn{K}) is used.
}
  \item{theta.k}{
A grid of theta values for the ability distribution.
}
  \item{est.a.item}{
Should item parameters \eqn{a_i} be estimated?
}
  \item{est.c.rater}{
Type of estimation for item-rater parameters \eqn{c_{ir}}
in the signal detection model. Options are \code{'n'} (no estimation),
\code{'e'} (set all parameters equal to each other),
\code{'i'} (itemwise estimation), \code{'r'} (rater wise estimation)
and \code{'a'} (all parameters are estimated independently
from each other).
}
  \item{est.d.rater}{
Type of estimation of \eqn{d} parameters. Options are
the same as in \code{est.c.rater}.
}
\item{est.mean}{Optional logical indicating whether the mean of the
trait distribution should be estimated.}
\item{est.sigma}{Optional logical indicating whether the standard deviation of the
trait distribution should be estimated.}
\item{skillspace}{Specified \eqn{\theta} distribution type. It can be
    \code{"normal"} or \code{"discrete"}. In the latter case,
    all probabilities of the distribution are separately
       estimated.}
\item{tau.item.fixed}{Optional matrix with three columns specifying
fixed \eqn{\tau} parameters. The first two columns denote item and
category indices, the third the fixed value. See Example 3.}
\item{a.item.fixed}{Optional matrix with two columns specifying fixed
    \eqn{a} parameters. First column: Item index. Second column:
Fixed \eqn{a} parameter.}
  \item{d.min}{
Minimal \eqn{d} parameter to be estimated
}
  \item{d.max}{
Maximal \eqn{d} parameter to be estimated
}
\item{d.start}{
Starting value(s) of \eqn{d} parameters
}
\item{c.start}{
Starting values of \eqn{c} parameters
}
\item{tau.start}{
Starting values of \eqn{\tau} parameters
}
\item{sd.start}{
Starting value for trait standard deviation
}
\item{d.prior}{Normal prior \eqn{N(M,S^2)} for \eqn{d} parameters}
\item{c.prior}{Normal prior for \eqn{c} parameters. The prior for
    parameter \eqn{c_{irk}} is defined as \eqn{M \cdot ( k - 0.5) }
    where \eqn{M} is \code{c.prior[1]}. }
\item{tau.prior}{Normal prior for \eqn{\tau} parameters}
\item{a.prior}{Normal prior for \eqn{a} parameters}
\item{link_item}{Type of item response function for latent responses.
Can be \code{"GPCM"} for the generalized partial credit model or
\code{"GRM"} for the graded response model.
}
  \item{max.increment}{
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
\item{optimizer}{Choice of optimization function in M-step for
item parameters. Options are \code{"nlminb"} for \code{\link[stats:nlminb]{stats::nlminb}}
and \code{"optim"} for \code{\link[stats:optim]{stats::optim}}.
}
  \item{object}{
Object of class \code{rm.sdt}
}
\item{file}{Optional file name in which summary should be written.}
\item{x}{Object of class \code{rm.sdt}}
\item{ask}{Optional logical indicating whether a new plot should be asked for.}
\item{type}{Factor score estimation method. Up to now,
    only \code{type="EAP"} is supported.}
\item{\dots}{
Further arguments to be passed
}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETAILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\details{
The specification of the model follows DeCarlo et al. (2011).
The second level models the ideal rating (latent response) \eqn{\eta=0, ...,K}
of person \eqn{p} on item \eqn{i}. The option \code{link_item='GPCM'} follows the
generalized partial credit model
\deqn{ P( \eta_{pi}=\eta | \theta_p ) \propto
exp( a_{i} q_{i \eta } \theta_p - \tau_{i \eta } ) }. The option \code{link_item='GRM'} employs the
graded response model   \deqn{ P( \eta_{pi}=\eta | \theta_p )=
\Psi( \tau_{i,\eta + 1} - a_i \theta_p ) - \Psi( \tau_{i,\eta} - a_i \theta_p ) }

At the first level, the ratings \eqn{X_{pir}} for
person \eqn{p} on item \eqn{i} and rater \eqn{r}
are modeled as a signal detection model
\deqn{ P( X_{pir} \le k | \eta_{pi} )=
G( c_{irk} - d_{ir} \eta_{pi} )}
where \eqn{G} is the logistic distribution function
and the categories are \eqn{k=1,\ldots, K+1}.
Note that the item response model can be equivalently written
as \deqn{ P( X_{pir} \ge k | \eta_{pi} )=
G(   d_{ir} \eta_{pi} - c_{irk})}

The thresholds \eqn{c_{irk}} can be further restricted to
\eqn{c_{irk}=c_{k}} (\code{est.c.rater='e'}),
\eqn{c_{irk}=c_{ik}} (\code{est.c.rater='i'}) or
\eqn{c_{irk}=c_{ir}} (\code{est.c.rater='r'}). The same
holds for rater precision parameters \eqn{d_{ir}}.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\value{
A list with following entries:
\item{deviance}{Deviance}
\item{ic}{Information criteria and number of parameters}
\item{item}{Data frame with item parameters. The columns
\code{N} and \code{M} denote the number of observed ratings and the
observed mean of all ratings, respectively. \cr
In addition to item parameters \eqn{\tau_{ik}} and \eqn{a_i}, the mean
for the latent response (\code{latM}) is computed as
\eqn{E( \eta_i )=\sum_p P( \theta_p ) q_{ik} P( \eta_i=k | \theta_p ) }
which provides an item parameter at the original metric of ratings. The latent standard
deviation (\code{latSD}) is computed in the same manner.
}
\item{rater}{Data frame with rater parameters.
Transformed \eqn{c} parameters
(\code{c_x.trans}) are computed as \eqn{c_{irk} / ( d_{ir}  )}.
}
\item{person}{Data frame with person parameters: EAP and corresponding standard
errors}
\item{EAP.rel}{EAP reliability}
\item{EAP.rel}{EAP reliability}
\item{mu}{Mean of the trait distribution}
\item{sigma}{Standard deviation of the trait distribution}
\item{tau.item}{Item parameters \eqn{\tau_{ik}}}
\item{se.tau.item}{Standard error of item parameters \eqn{\tau_{ik}}}
\item{a.item}{Item slopes \eqn{a_i}}
\item{se.a.item}{Standard error of item slopes \eqn{a_i}}
\item{c.rater}{Rater parameters \eqn{c_{irk}}}
\item{se.c.rater}{Standard error of rater severity parameter \eqn{c_{irk}}}
\item{d.rater}{Rater slope parameter \eqn{d_{ir}}}
\item{se.d.rater}{Standard error of rater slope parameter \eqn{d_{ir}}}
\item{f.yi.qk}{Individual likelihood}
\item{f.qk.yi}{Individual posterior distribution}
\item{probs}{Item probabilities at grid \code{theta.k}. Note that these
  probabilities are calculated on the pseudo items \eqn{i \times r},
  i.e. the interaction of item and rater.}
\item{prob.item}{Probabilities \eqn{P( \eta_i=\eta | \theta )} of
    latent item responses evaluated at theta grid \eqn{\theta_p}.
}
\item{n.ik}{Expected counts}
\item{pi.k}{Estimated trait distribution \eqn{P(\theta_p)}.
}
\item{maxK}{Maximum number of categories}
\item{procdata}{Processed data}
\item{iter}{Number of iterations}
\item{\dots}{Further values}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\references{
DeCarlo, L. T. (2005). A model of rater behavior in essay grading
based on signal detection theory.
\emph{Journal of Educational Measurement}, \bold{42}, 53-76.

DeCarlo, L. T. (2010). \emph{Studies of a latent-class signal-detection model
for constructed response scoring II: Incomplete and hierarchical designs}.
ETS Research Report ETS RR-10-08. Princeton NJ: ETS.

DeCarlo, T., Kim, Y., & Johnson, M. S. (2011).
A hierarchical rater model for constructed responses,
with a signal detection rater model.
\emph{Journal of Educational Measurement},
\bold{48}, 333-356.

Robitzsch, A., & Steinfeld, J. (2018). Item response models for human ratings: Overview,
estimation methods, and implementation in R.
\emph{Psychological Test and Assessment Modeling, 60}(1), 101-139.

Vermunt, J. K. (2008). Latent class and finite mixture models for
multilevel data sets. \emph{Statistical Methods in Medical Research, 17}, 33-51.
}


\author{
Alexander Robitzsch
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
The facets rater model can be estimated with \code{\link{rm.facets}}.
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\examples{
#############################################################################
# EXAMPLE 1: Hierarchical rater model (HRM-SDT) data.ratings1
#############################################################################
data(data.ratings1)
dat <- data.ratings1

\dontrun{
# Model 1: Partial Credit Model: no rater effects
mod1 <- sirt::rm.sdt( dat[, paste0( "k",1:5) ], rater=dat$rater,
            pid=dat$idstud, est.c.rater="n", d.start=100,  est.d.rater="n" )
summary(mod1)

# Model 2: Generalized Partial Credit Model: no rater effects
mod2 <- sirt::rm.sdt( dat[, paste0( "k",1:5) ], rater=dat$rater,
            pid=dat$idstud, est.c.rater="n", est.d.rater="n",
            est.a.item=TRUE, d.start=100)
summary(mod2)

# Model 3: Equal effects in SDT
mod3 <- sirt::rm.sdt( dat[, paste0( "k",1:5) ], rater=dat$rater,
            pid=dat$idstud, est.c.rater="e", est.d.rater="e")
summary(mod3)

# Model 4: Rater effects in SDT
mod4 <- sirt::rm.sdt( dat[, paste0( "k",1:5) ], rater=dat$rater,
            pid=dat$idstud, est.c.rater="r", est.d.rater="r")
summary(mod4)

#############################################################################
# EXAMPLE 2: HRM-SDT data.ratings3
#############################################################################

data(data.ratings3)
dat <- data.ratings3
dat <- dat[ dat$rater < 814, ]
psych::describe(dat)

# Model 1: item- and rater-specific effects
mod1 <- sirt::rm.sdt( dat[, paste0( "crit",c(2:4)) ], rater=dat$rater,
            pid=dat$idstud, est.c.rater="a", est.d.rater="a" )
summary(mod1)
plot(mod1)

# Model 2: Differing number of categories per variable
mod2 <- sirt::rm.sdt( dat[, paste0( "crit",c(2:4,6)) ], rater=dat$rater,
            pid=dat$idstud, est.c.rater="a", est.d.rater="a")
summary(mod2)
plot(mod2)

#############################################################################
# EXAMPLE 3: Hierarchical rater model with discrete skill spaces
#############################################################################

data(data.ratings3)
dat <- data.ratings3
dat <- dat[ dat$rater < 814, ]
psych::describe(dat)

# Model 1: Discrete theta skill space with values of 0,1,2 and 3
mod1 <- sirt::rm.sdt( dat[, paste0( "crit",c(2:4)) ], theta.k=0:3, rater=dat$rater,
            pid=dat$idstud, est.c.rater="a", est.d.rater="a", skillspace="discrete" )
summary(mod1)
plot(mod1)

# Model 2: Modelling of one item by using a discrete skill space and
#          fixed item parameters

# fixed tau and a parameters
tau.item.fixed <- cbind( 1, 1:3,  100*cumsum( c( 0.5, 1.5, 2.5)) )
a.item.fixed <- cbind( 1, 100 )
# fit HRM-SDT
mod2 <- sirt::rm.sdt( dat[, "crit2", drop=FALSE], theta.k=0:3, rater=dat$rater,
            tau.item.fixed=tau.item.fixed,a.item.fixed=a.item.fixed, pid=dat$idstud,
            est.c.rater="a", est.d.rater="a", skillspace="discrete" )
summary(mod2)
plot(mod2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Rater model}
\keyword{Signal detection model}
\keyword{summary}
\keyword{plot}
