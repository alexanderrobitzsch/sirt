%% File Name: rasch.mml.Rd
%% File Version: 4.006

\name{rasch.mml2}
\alias{rasch.mml2}
\alias{summary.rasch.mml}
\alias{plot.rasch.mml}
\alias{logLik.rasch.mml}
\alias{anova.rasch.mml}
\alias{IRT.irfprob.rasch.mml}
\alias{IRT.likelihood.rasch.mml}
\alias{IRT.posterior.rasch.mml}
\alias{IRT.modelfit.rasch.mml}
\alias{IRT.expectedCounts.rasch.mml}
\alias{summary.IRT.modelfit.rasch.mml}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of the Generalized Logistic Item Response Model,
Ramsay's Quotient Model, Nonparametric Item Response Model,
Pseudo-Likelihood Estimation and a Missing Data Item Response Model
}
\description{
This function employs marginal maximum likelihood estimation
of item response models for dichotomous data.
First, the Rasch type model (generalized
item response model) can be estimated. The generalized logistic
link function (Stukel, 1988) can be estimated or fixed for conducting
IRT with different link functions than the logistic one. The Four-Parameter
logistic item response model is a special case of this model
(Loken & Rulison, 2010). Second, Ramsay's quotient model (Ramsay, 1989)
can be estimated by specifying \code{irtmodel="ramsay.qm"}.
Third, quite general item response functions can be estimated
in a nonparametric framework (Rossi, Wang & Ramsay, 2002).
Fourth, pseudo-likelihood estimation for fractional item responses can be
conducted for Rasch type models. Fifth, a simple two-dimensional
missing data item response model (\code{irtmodel='missing1'};
Mislevy & Wu, 1996) can be estimated.

See Details for more explanations.
}
\usage{
rasch.mml2( dat, theta.k=seq(-6,6,len=21), group=NULL, weights=NULL,
   constraints=NULL, glob.conv=10^(-5), parm.conv=10^(-4), mitermax=4,
   mmliter=1000, progress=TRUE,  fixed.a=rep(1,ncol(dat)),
   fixed.c=rep(0,ncol(dat)), fixed.d=rep(1,ncol(dat)),
   fixed.K=rep(3,ncol(dat)), b.init=NULL, est.a=NULL, est.b=NULL,
   est.c=NULL, est.d=NULL, min.b=-99, max.b=99, min.a=-99, max.a=99,
   min.c=0, max.c=1, min.d=0, max.d=1, prior.b=NULL, prior.a=NULL, prior.c=NULL,
   prior.d=NULL, est.K=NULL, min.K=1, max.K=20, beta.init=NULL, min.beta=-8,
   pid=1:(nrow(dat)), trait.weights=NULL,  center.trait=TRUE, center.b=FALSE, alpha1=0,
   alpha2=0,est.alpha=FALSE, equal.alpha=FALSE, designmatrix=NULL, alpha.conv=parm.conv,
   numdiff.parm=0.00001, numdiff.alpha.parm=numdiff.parm, distribution.trait="normal",
   Qmatrix=NULL, variance.fixed=NULL, variance.init=NULL,
   mu.fixed=cbind(seq(1,ncol(Qmatrix)),rep(0,ncol(Qmatrix))),
   irtmodel="raschtype", npformula=NULL, npirt.monotone=TRUE,
   use.freqpatt=is.null(group), delta.miss=0, est.delta=rep(NA,ncol(dat)),
   ... )

\method{summary}{rasch.mml}(object, file=NULL, ...)

\method{plot}{rasch.mml}(x, items=NULL, xlim=NULL, main=NULL, ...)

\method{anova}{rasch.mml}(object,...)

\method{logLik}{rasch.mml}(object,...)

\method{IRT.irfprob}{rasch.mml}(object,...)

\method{IRT.likelihood}{rasch.mml}(object,...)

\method{IRT.posterior}{rasch.mml}(object,...)

\method{IRT.modelfit}{rasch.mml}(object,...)

\method{IRT.expectedCounts}{rasch.mml}(object,...)

\method{summary}{IRT.modelfit.rasch.mml}(object,...)
}
%- maybe also 'usage' for other objects documented here.

%% Undocumented arguments in documentation object 'rasch.mml2'
%%  'Qmatrix' 'variance.fixed' 'mu.fixed'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\arguments{
  \item{dat}{
An \eqn{N \times I} data frame of dichotomous item responses. \cr
For the missing data item response model (\code{irtmodel='missing1'}),
code item responses by \code{9} which should be treated by the missing
  data model. Other missing responses can be coded by \code{NA}.
}
\item{theta.k}{
Optional vector of discretized theta values. For multidimensional
IRT models with \eqn{D} dimensions, it is a matrix with \eqn{D} columns.
}
  \item{group}{
Vector of integers with group identifiers in multiple group estimation.
The multiple group does not work for \code{irtmodel="missing1"}.
}
  \item{weights}{
Optional vector of person weights (sample weights).
}
  \item{constraints}{
Constraints on \code{b} parameters (item difficulties). It must be
a matrix with two columns: the first column contains item names, the
second column fixed parameter values.
}
  \item{glob.conv}{
Convergence criterion for deviance
}
  \item{parm.conv}{
Convergence criterion for item parameters
}
  \item{mitermax}{
Maximum number of iterations in M step. This argument does only
apply for the estimation of the \eqn{b} parameters.
}
  \item{mmliter}{
Maximum number of iterations
}
  \item{progress}{
Should progress be displayed at the console?
}
  \item{fixed.a}{
Fixed or initial \eqn{a} parameters
}
  \item{fixed.c}{
Fixed or initial \eqn{c} parameters
}
  \item{fixed.d}{
Fixed or initial \eqn{d} parameters
}
  \item{fixed.K}{
Fixed or initial \eqn{K} parameters in Ramsay's
quotient model.
}
  \item{b.init}{
Initial \eqn{b} parameters
}
  \item{est.a}{
Vector of integers which indicate which \eqn{a}
parameters should be estimated. Equal integers correspond
to the same estimated parameters.
}
  \item{est.b}{
Vector of integers which indicate which \eqn{b}
parameters should be estimated. Equal integers correspond
to the same estimated parameters.
}
  \item{est.c}{
Vector of integers which indicate which \eqn{c}
parameters should be estimated. Equal integers correspond
to the same estimated parameters.
}
  \item{est.d}{
Vector of integers which indicate which \eqn{d}
parameters should be estimated. Equal integers correspond
to the same estimated parameters.
}
  \item{min.b}{
Minimal \eqn{b} parameter to be estimated
}
  \item{max.b}{
Maximal \eqn{b} parameter to be estimated
}
  \item{min.a}{
Minimal \eqn{a} parameter to be estimated
}
  \item{max.a}{
Maximal \eqn{a} parameter to be estimated
}
  \item{min.c}{
Minimal \eqn{c} parameter to be estimated
}
  \item{max.c}{
Maximal \eqn{c} parameter to be estimated
}
  \item{min.d}{
Minimal \eqn{d} parameter to be estimated
}
  \item{max.d}{
Maximal \eqn{d} parameter to be estimated
}
\item{prior.b}{Optional prior distribution for \eqn{b} parameters:
\eqn{N(\mu, \sigma)}. Input is a vector of length two with parameters
\eqn{\mu} and \eqn{\sigma}.}
\item{prior.a}{Optional prior distribution for \eqn{a} parameters:
\eqn{N(\mu, \sigma)}. Input is a vector of length two with parameters
\eqn{\mu} and \eqn{\sigma}.}
\item{prior.c}{Optional prior distribution for \eqn{c} parameters:
\eqn{Beta(a, b)}. Input is a vector of length two with parameters
\eqn{a} and \eqn{b}.}
\item{prior.d}{Optional prior distribution for \eqn{d} parameters:
\eqn{Beta(a, b)}. Input is a vector of length two with parameters
\eqn{a} and \eqn{b}.}
  \item{est.K}{
Vector of integers which indicate which \eqn{K}
parameters should be estimated. Equal integers correspond
to the same estimated parameters.
}
  \item{min.K}{
Minimal \eqn{K} parameter to be estimated
}
  \item{max.K}{
Maximal \eqn{K} parameter to be estimated
}
\item{beta.init}{Optional vector of initial \eqn{\beta} parameters}
\item{min.beta}{Minimum \eqn{\beta} parameter to be estimated.}
  \item{pid}{
Optional vector of person identifiers
}
  \item{trait.weights}{
Optional vector of trait weights for a fixing
the trait distribution.
}
  \item{center.trait}{
Should the trait distribution be centered
}
\item{center.b}{An optional logical indicating whether \eqn{b} parameters
should be centered at each dimension}
  \item{alpha1}{
Fixed or initial \eqn{\alpha_1} parameter
}
  \item{alpha2}{
Fixed or initial \eqn{\alpha_2} parameter
}
  \item{est.alpha}{
Should \eqn{\alpha} parameters be estimated?
}
  \item{equal.alpha}{
Estimate \eqn{\alpha} parameters under the
assumption \eqn{\alpha_1=\alpha_2}?
}
  \item{designmatrix}{
Design matrix for item difficulties \eqn{b} to estimate
linear logistic test models
}
  \item{alpha.conv}{
Convergence criterion for \eqn{\alpha} parameter
}
  \item{numdiff.parm}{
Parameter for numerical differentiation
}
  \item{numdiff.alpha.parm}{
Parameter for numerical differentiation for \eqn{\alpha}
parameter
}
  \item{distribution.trait}{
Assumed trait distribution. The default is the normal
distribution (\code{"normal"}). Log-linear smoothing of the
trait distribution is also possible (\code{"smooth2"},
\code{"smooth3"} or \code{"smooth4"} for smoothing up to
2, 3 or 4 moments, respectively).
}
  \item{Qmatrix}{
The Q-matrix
}
  \item{variance.fixed}{
Matrix for fixing covariance matrix (See Examples)
}
\item{variance.init}{Optional initial covariance matrix}
  \item{mu.fixed}{
Matrix for fixing mean vector (See Examples)
}

  \item{irtmodel}{
Specify estimable IRT models: \code{raschtype} (Rasch type model),
\code{ramsay.qm} (Ramsay's quotient model), \code{npirt} (Nonparametric
item response model). If \code{npirt} is used as the argument
for \code{irtmodel}, the argument \code{npformula}
specifies different item response functions in the
\R formula framework (like \code{"y~I(theta^2)"}; see Examples).
For estimating the missing data item response model, use
\code{irtmodel='missing1'}.
}
  \item{npformula}{
A string or a vector which contains \R formula objects for specifying
the item response function. For example, \code{"y~theta"} is the specification
of the 2PL model (see Details). If \code{irtmodel="npirt"} and \code{npformula}
is not specified, then an unrestricted item response functions on the
grid of \eqn{\theta} values is estimated.
}

\item{npirt.monotone}{Should nonparametrically estimated item response functions
    be monotone? The default is \code{TRUE}. This function applies only
    to \code{irtmodel='npirt'} and \code{npformula=NULL}.
}
\item{use.freqpatt}{A logical if frequencies of pattern should be used or not.
    The default is \code{is.null(group)}. This means that for single
    group analyses, frequency patterns are used but not for multiple
    groups. If data processing times are large, then \code{use.freqpatt=FALSE}
    is recommended.}
\item{delta.miss}{Missingness parameter \eqn{\delta} quantifying the meaning
    of responding to an item between the two extremes of ignoring
     missing responses and setting all missing responses to incorrect}
\item{est.delta}{Vector with indices indicating the \eqn{\delta} parameters
    to be estimated if \code{irtmodel="missing1"}.}
  \item{object}{
Object of class \code{rasch.mml}
}
\item{x}{Object of class \code{rasch.mml}}
\item{items}{Vector of integer or item names which should be plotted}
\item{xlim}{Specification for \code{xlim} in plot}
\item{main}{Title of the plot}
\item{file}{Optional file name for summary output}
  \item{\dots}{
Further arguments to be passed
}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DETAILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\details{
The item response function of the generalized item response model
(\code{irtmodel="raschtype"}; Stukel, 1988) can be written as
\deqn{P( X_{pi}=1 | \theta_{pd} )=c_i + (d_i - c_i ) g_{\alpha_1, \alpha_2}
[ a_i ( \theta_{pd} - b_i ) ] }
where \eqn{g} is the generalized logistic link function depending
on parameters \eqn{\alpha_1} and \eqn{\alpha_2}.

For the most important link functions the specifications are (Stukel, 1988):

logistic link function: \eqn{\alpha_1=0} and \eqn{\alpha_2=0} \cr
probit link function: \eqn{\alpha_1=0.165} and \eqn{\alpha_2=0.165} \cr
loglog link function: \eqn{\alpha_1=-0.037} and \eqn{\alpha_2=0.62} \cr
cloglog link function: \eqn{\alpha_1=0.62} and \eqn{\alpha_2=-0.037}

See \code{\link{pgenlogis}} for exact transformation formulas of
the mentioned link functions. \cr

A \eqn{D}-dimensional model can also be specified
but only allows for between item dimensionality
(one item loads on one and only dimension).
Setting \eqn{c_i=0}, \eqn{d_i=1} and  \eqn{a_i=1} for all items \eqn{i},
an additive item response model
\deqn{P( X_{pi}=1 | \theta_p )=g_{\alpha_1, \alpha_2} ( \theta_p - b_i  ) }
is estimated.

Ramsay's quotient model (\code{irtmodel="qm.ramsay"}) uses
the item response function
\deqn{P( X_{pi}=1 | \theta_p )=\frac{ \exp(\theta_p / b_i)}
    { K_i + \exp (\theta_p / b_i )} }

Quite general unidimensional item response models can be estimated
in a nonparametric framework (\code{irtmodel="npirt"}). The response
functions are a linear combination of transformed \eqn{\theta}
values
\deqn{logit[ P( X_{pi}=1 | \theta_p ) ]=Y_\theta \beta }
Where \eqn{Y_\theta} is a design matrix of \eqn{\theta} and
\eqn{\beta} are item parameters to be estimated.
The formula \eqn{Y_\theta \beta} can be specified in the \R formula
framework (see Example 3, Model 3c).

Pseudo-likelihood estimation can be conducted for fractional item response data
as input (i.e. some item response \eqn{x_{pi}} do have values
between 0 and 1). Then the pseudo-likelihood \eqn{L_p} for person \eqn{p} is defined as
\deqn{ L_p=\prod_i P_i ( \theta_p )^{x_{pi}} [1-P_i ( \theta_p )]^{(1-x_{pi})}}
Note that for dichotomous responses this term corresponds to the ordinary
likelihood. See Example 7.

A special two-dimensional missing data item response model (\code{irtmodel="missing1"})
is implemented according to Mislevy and Wu (1996).
Besides an unidimensional ability \eqn{\theta_p},
an individual response propensity \eqn{\xi_p} is proposed. We define
item responses \eqn{X_{pi}} and response indicators \eqn{R_{pi}} indicating whether
item responses \eqn{X_{pi}} are observed or not. Denoting the logistic function
by \eqn{L}, the item response model for ability is defined as
\deqn{ P( X_{pi}=1  | \theta_p, \xi_p )=P( X_{pi}=1 | \theta_p )
=L( \theta_p - b_i )}
We also define a measurement model for response indicators \eqn{R_{pi}} which depends
on the item response \eqn{X_{pi}} itself:
\deqn{P( R_{pi}=1 | X_{pi}=k, \theta_p, \xi_p )=
    P( R_{pi}=1 | X_{pi}=k, \xi_p )=
L \left[  \xi_p -  \beta_i - k \delta _i  \right] \quad \mbox{ for }
\quad k=0,1}
If \eqn{\delta _i=0}, then the probability of responding to an item is independent
of the incompletely observed item \eqn{X_{pi}} which is an
item response model with nonignorable missings (Holman & Glas, 2005;
see also Pohl, Graefe & Rose, 2014).
If \eqn{\delta _i} is a large negative number (e.g. \eqn{\delta=-100}), then
it follows \eqn{P( R_{pi}=1 | X_{pi}=1, \theta_p, \xi_p )=1}
and as a consequence it holds that \eqn{P(X_{pi}=1 | R_{pi}=0, \theta_p, \xi_p)=0},
which is equivalent to treating
all missing item responses as incorrect. The missingness parameter
\eqn{\delta} can be specified
by the user and studied as a sensitivity analysis under different
missing not at random assumptions or can be estimated by choosing
\code{est.delta=TRUE}.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\value{
A list with following entries
\item{dat}{Original data frame}
\item{item}{Estimated item parameters in the generalized
    item response model}
\item{item2}{Estimated item parameters for Ramsay's quotient model}
\item{trait.distr}{Discretized ability distribution points and probabilities}
\item{mean.trait}{Estimated mean vector}
\item{sd.trait}{Estimated standard deviations}
\item{skewness.trait}{Estimated skewnesses}
\item{deviance}{Deviance}
\item{pjk}{Estimated probabilities of item correct evaluated at \code{theta.k}}
\item{rprobs}{Item response probabilities like in \code{pjk}, but slightly
extended to accommodate all categories}
\item{person}{Person parameter estimates: mode (\code{MAP}) and
    mean (\code{EAP}) of the posterior distribution}
\item{pid}{Person identifier}
\item{ability.est.pattern}{Response pattern estimates}
\item{f.qk.yi}{Individual posterior distribution}
\item{f.yi.qk}{Individual likelihood}
\item{fixed.a}{Estimated \eqn{a} parameters}
\item{fixed.c}{Estimated \eqn{c} parameters}
\item{G}{Number of groups}
\item{alpha1}{Estimated \eqn{\alpha_1} parameter in generalized logistic
    item response model}
\item{alpha2}{Estimated \eqn{\alpha_2} parameter in generalized logistic
    item response model}
\item{se.b}{Standard error of \eqn{b} parameter in generalized logistic model
    or Ramsay's quotient model}
\item{se.a}{Standard error of \eqn{a} parameter in generalized logistic model}
\item{se.c}{Standard error of \eqn{c} parameter in generalized logistic model}
\item{se.d}{Standard error of \eqn{d} parameter in generalized logistic model}
\item{se.alpha}{Standard error of \eqn{\alpha} parameter in generalized
    logistic model}
\item{se.K}{Standard error of \eqn{K} parameter in Ramsay's quotient model}
\item{iter}{Number of iterations}
\item{reliability}{EAP reliability}
\item{irtmodel}{Type of estimated item response model}
\item{D}{Number of dimensions}
\item{mu}{Mean vector (for multidimensional models)}
\item{Sigma.cov}{Covariance matrix (for multdimensional models)}
\item{theta.k}{Grid of discretized ability distributions}
\item{trait.weights}{Fixed vector of probabilities for the ability distribution}
\item{pi.k}{Trait distribution}
\item{ic}{Information criteria}
\item{esttype}{Estimation type: \code{ll} (Log-Likelihood),
    \code{pseudoll} (Pseudo-Log-Likelihood)
}
\item{\dots}{}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  REFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\references{
Holman, R., & Glas, C. A. (2005). Modelling non-ignorable missing-data mechanisms
with item response theory models.
\emph{British Journal of Mathematical and Statistical Psychology, 58}(1), 1-17.
\doi{10.1348/000711005X47168}

Loken, E., & Rulison, K. L. (2010). Estimation of a four-parameter
item response theory model. \emph{British Journal of Mathematical
and Statistical Psychology, 63}(3), 509-525.
\doi{10.1348/000711009X474502}

Mislevy, R. J., & Wu, P. K. (1996). \emph{Missing responses and IRT ability
estimation: Omits, choice, time Limits, and adaptive testing}.
ETS Research Report ETS RR-96-30. Princeton, ETS.
\doi{10.1002/j.2333-8504.1996.tb01708.x}

Pohl, S., Graefe, L., & Rose, N. (2014). Dealing with omitted and
not-reached items in competence tests evaluating approaches accounting for
missing responses in item response theory models.
\emph{Educational and Psychological Measurement, 74}(3), 423-452.
\doi{10.1177/0013164413504926}

Ramsay, J. O. (1989). A comparison of three simple test theory models.
\emph{Psychometrika, 54}, 487-499.
\doi{10.1007/BF02294631}

Rossi, N., Wang, X., & Ramsay, J. O. (2002). Nonparametric item response
function estimates with the EM algorithm.
\emph{Journal of Educational and Behavioral Statistics, 27}(3), 291-317.
\doi{10.3102/10769986027003291}

Stukel, T. A. (1988). Generalized logistic models.
\emph{Journal of the American Statistical Association, 83}(402), 426-431.
\doi{10.1080/01621459.1988.10478613}

van der Maas, H. J. L., Molenaar, D., Maris, G., Kievit, R. A., &
Borsboom, D. (2011).
Cognitive psychology meets psychometric theory: On the relation between
process models for decision making and latent variable models for
individual differences.
\emph{Psychological Review, 118}(2), 339-356.
\doi{10.1037/a0022749}
}


\note{
Multiple group estimation is not possible for
Ramsay's quotient model and multdimensional models.
}

%Be careful with this function, especially with multiple groups.
%It is just in an experimental stage.

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Simulate the generalized logistic Rasch model with \code{\link{sim.raschtype}}.

Simulate Ramsay's quotient model with \code{\link{sim.qm.ramsay}}.

Simulate locally dependent item response data using \code{\link{sim.rasch.dep}}.

For an assessment of global model fit see \code{\link{modelfit.sirt}}.

See \code{\link[CDM:itemfit.sx2]{CDM::itemfit.sx2}} for item fit
statistics.
}

% \code{\link[mice:mice]{mice}} (\pkg{mice})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



\examples{
#############################################################################
# EXAMPLE 1: Reading dataset
#############################################################################

library(CDM)
data(data.read)
dat <- data.read
I <- ncol(dat) # number of items

# Rasch model
mod1 <- sirt::rasch.mml2( dat )
summary(mod1)
plot( mod1 )    # plot all items
# title 'Rasch model', display curves from -3 to 3 only for items 1, 5 and 8
plot(mod1, main="Rasch model Items 1, 5 and 8", xlim=c(-3,3), items=c(1,5,8) )

# Rasch model with constraints on item difficulties
# set item parameters of A1 and C3 equal to -2
constraints <- data.frame( c("A1","C3"), c(-2,-2) )
mod1a <- sirt::rasch.mml2( dat, constraints=constraints)
summary(mod1a)

# estimate equal item parameters for 1st and 11th item
est.b <- 1:I
est.b[11] <- 1
mod1b <- sirt::rasch.mml2( dat, est.b=est.b )
summary(mod1b)

# estimate Rasch model with skew trait distribution
mod1c <- sirt::rasch.mml2( dat, distribution.trait="smooth3")
summary(mod1c)

# 2PL model
mod2 <- sirt::rasch.mml2( dat, est.a=1:I )
summary(mod2)
plot(mod2)    # plot 2PL item response curves

# extract individual likelihood
llmod2 <- IRT.likelihood(mod2)
str(llmod2)

\dontrun{
library(CDM)
# model comparisons
CDM::IRT.compareModels(mod1, mod1c, mod2 )
anova(mod1,mod2)

# assess model fit
smod1 <- IRT.modelfit(mod1)
smod2 <- IRT.modelfit(mod2)
IRT.compareModels(smod1, smod2)

# set some bounds for a and b parameters
mod2a <- sirt::rasch.mml2( dat, est.a=1:I, min.a=.7, max.a=2, min.b=-2 )
summary(mod2a)

# 3PL model
mod3 <- sirt::rasch.mml2( dat, est.a=1:I, est.c=1:I,
              mmliter=400 # maximal 400 iterations
                 )
summary(mod3)

# 3PL model with fixed guessing paramters of .25 and equal slopes
mod4 <- sirt::rasch.mml2( dat, fixed.c=rep( .25, I )   )
summary(mod4)

# 3PL model with equal guessing paramters for all items
mod5 <- sirt::rasch.mml2( dat, est.c=rep(1, I )   )
summary(mod5)

# difficulty + guessing model
mod6 <- sirt::rasch.mml2( dat, est.c=1:I   )
summary(mod6)

# 4PL model
mod7 <- sirt::rasch.mml2( dat, est.a=1:I, est.c=1:I, est.d=1:I,
            min.d=.95, max.c=.25)
        # set minimal d and maximal c parameter to .95 and .25
summary(mod7)

# 4PL model with prior distributions
mod7b <- sirt::rasch.mml2( dat, est.a=1:I, est.c=1:I, est.d=1:I, prior.a=c(1,2),
            prior.c=c(5,17), prior.d=c(20,2) )
summary(mod7b)

# constrained 4PL model
# equal slope, guessing and slipping parameters
mod8 <- sirt::rasch.mml2( dat,est.c=rep(1,I), est.d=rep(1,I) )
summary(mod8)

# estimation of an item response model with an
# uniform theta distribution
theta.k <- seq( 0.01, .99, len=20 )
trait.weights <- rep( 1/length(theta.k), length(theta.k) )
mod9 <- sirt::rasch.mml2( dat, theta.k=theta.k, trait.weights=trait.weights,
              normal.trait=FALSE, est.a=1:12  )
summary(mod9)

#############################################################################
# EXAMPLE 2: Longitudinal data
#############################################################################

data(data.long)
dat <- data.long[,-1]

# define Q loading matrix
Qmatrix <- matrix( 0, 12, 2 )
Qmatrix[1:6,1] <- 1 # T1 items
Qmatrix[7:12,2] <- 1    # T2 items

# define restrictions on item difficulties
est.b <- c(1,2,3,4,5,6,   3,4,5,6,7,8)
mu.fixed <- cbind(1,0)
    # set first mean to 0 for identification reasons

# Model 1: 2-dimensional Rasch model
mod1 <- sirt::rasch.mml2( dat, Qmatrix=Qmatrix, miterstep=4,
            est.b=est.b,  mu.fixed=mu.fixed, mmliter=30 )
summary(mod1)
plot(mod1)
##     Plot function is only applicable for unidimensional models
}

#############################################################################
# EXAMPLE 3: One group, estimation of alpha parameter in the generalized logistic model
#############################################################################

# simulate theta values
set.seed(786)
N <- 1000                  # number of persons
theta <- stats::rnorm( N, sd=1.5 ) # N persons with SD 1.5
b <- seq( -2, 2, len=15)

# simulate data
dat <- sirt::sim.raschtype( theta=theta, b=b, alpha1=0, alpha2=-0.3 )

#  estimating alpha parameters
mod1 <- sirt::rasch.mml2( dat, est.alpha=TRUE, mmliter=30 )
summary(mod1)
plot(mod1)

\dontrun{
# fixed alpha parameters
mod1b <- sirt::rasch.mml2( dat, est.alpha=FALSE, alpha1=0, alpha2=-.3 )
summary(mod1b)

# estimation with equal alpha parameters
mod1c <- sirt::rasch.mml2( dat, est.alpha=TRUE, equal.alpha=TRUE )
summary(mod1c)

# Ramsay QM
mod2a <- sirt::rasch.mml2( dat, irtmodel="ramsay.qm" )
summary(mod2a)
}

# Ramsay QM with estimated K parameters
mod2b <- sirt::rasch.mml2( dat, irtmodel="ramsay.qm", est.K=1:15, mmliter=30)
summary(mod2b)
plot(mod2b)

\dontrun{
# nonparametric estimation of monotone item response curves
mod3a <- sirt::rasch.mml2( dat, irtmodel="npirt", mmliter=100,
            theta.k=seq( -3, 3, len=10) ) # evaluations at 10 theta grid points
# nonparametric ICC of first 4 items
round( t(mod3a$pjk)[1:4,], 3 )
summary(mod3a)
plot(mod3a)

# nonparametric IRT estimation without monotonicity assumption
mod3b <- sirt::rasch.mml2( dat, irtmodel="npirt", mmliter=10,
            theta.k=seq( -3, 3, len=10), npirt.monotone=FALSE)
plot(mod3b)

# B-Spline estimation of ICCs
library(splines)
mod3c <- sirt::rasch.mml2( dat, irtmodel="npirt",
             npformula="y~bs(theta,df=3)", theta.k=seq(-3,3,len=15) )
summary(mod3c)
round( t(mod3c$pjk)[1:6,], 3 )
plot(mod3c)

# estimation of quadratic item response functions: ~ theta + I( theta^2)
mod3d <- sirt::rasch.mml2( dat, irtmodel="npirt",
             npformula="y~theta + I(theta^2)" )
summary(mod3d)
plot(mod3d)

# estimation of a stepwise ICC function
# ICCs are constant on the theta domains: [-Inf,-1], [-1,1], [1,Inf]
mod3e <- sirt::rasch.mml2( dat, irtmodel="npirt",
             npformula="y~I(theta>-1 )+I(theta>1)" )
summary(mod3e)
plot(mod3e, xlim=c(-2.5,2.5) )

# 2PL model
mod4 <- sirt::rasch.mml2( dat,  est.a=1:15)
summary(mod4)

#############################################################################
# EXAMPLE 4: Two groups, estimation of generalized logistic model
#############################################################################

# simulate generalized logistic Rasch model in two groups
set.seed(8765)
N1 <- 1000     # N1=1000 persons in group 1
N2 <- 500      # N2=500 persons in group 2
dat1 <- sirt::sim.raschtype( theta=stats::rnorm( N1, sd=1.5 ), b=b,
            alpha1=-0.3, alpha2=0)
dat2 <- sirt::sim.raschtype( theta=stats::rnorm( N2, mean=-.5, sd=.75),
            b=b, alpha1=-0.3, alpha2=0)
dat1 <- rbind( dat1, dat2 )
group <- c( rep(1,N1), rep(2,N2))

mod1 <-  sirt::rasch.mml2( dat1, parm.conv=.0001, group=group, est.alpha=TRUE )
summary(mod1)

#############################################################################
# EXAMPLE 5: Multidimensional model
#############################################################################

#***
# (1) simulate data
set.seed(785)
library(mvtnorm)
N <- 500
theta <- mvtnorm::rmvnorm( N,mean=c(0,0), sigma=matrix( c(1.45,.5,.5,1.7), 2, 2 ))
I <- 10
# 10 items load on the first dimension
p1 <- stats::plogis( outer( theta[,1], seq( -2, 2, len=I ), "-" ) )
resp1 <- 1 * ( p1 > matrix( stats::runif( N*I ), nrow=N, ncol=I ) )
# 10 items load on the second dimension
p1 <- stats::plogis( outer( theta[,2], seq( -2, 2, len=I ), "-" ) )
resp2 <- 1 * ( p1 > matrix( stats::runif( N*I ), nrow=N, ncol=I ) )
#Combine the two sets of items into one response matrix
resp <- cbind(resp1,resp2)
colnames(resp) <- paste("I", 1:(2*I), sep="")
dat <- resp

# define Q-matrix
Qmatrix <- matrix( 0, 2*I, 2 )
Qmatrix[1:I,1] <- 1
Qmatrix[1:I+I,2] <- 1

#***
# (2) estimation of models
# 2-dimensional Rasch model
mod1 <- sirt::rasch.mml2( dat, Qmatrix=Qmatrix )
summary(mod1)

# 2-dimensional 2PL model
mod2 <- sirt::rasch.mml2( dat, Qmatrix=Qmatrix, est.a=1:(2*I) )
summary(mod2)

# estimation with some fixed variances and covariances
# set variance of 1st dimension to 1 and
#  covariance to zero
variance.fixed <- matrix( cbind(c(1,1), c(1,2), c(1,0)),
             byrow=FALSE, ncol=3 )
mod3 <- sirt::rasch.mml2( dat, Qmatrix=Qmatrix, variance.fixed=variance.fixed )
summary(mod3)

# constraints on item difficulties
#  useful for example in longitudinal linking
est.b <- c( 1:I, 1:I )
    # equal indices correspond to equally estimated item parameters
mu.fixed <- cbind( 1, 0 )
mod4 <- sirt::rasch.mml2( dat, Qmatrix=Qmatrix, est.b=est.b, mu.fixed=mu.fixed )
summary(mod4)

#############################################################################
# EXAMPLE 6: Two booklets with same items but with item context effects.
# Therefore, item slopes and item difficulties are assumed to be shifted in the
# second design group.
#############################################################################

#***
# simulate data
set.seed(987)
I <- 10     # number of items
# define person design groups 1 and 2
n1 <- 700
n2 <- 1500
# item difficulties group 1
b1 <- seq(-1.5,1.5,length=I)
# item slopes group 1
a1 <- rep(1, I)
# simulate data group 1
dat1 <- sirt::sim.raschtype( stats::rnorm(n1), b=b1, fixed.a=a1 )
colnames(dat1) <- paste0("I", 1:I, "des1" )
# group 2
b2 <- b1 - .15
a2 <- 1.1*a1
# Item parameters are slightly transformed in the second group
# compared to the first group. This indicates possible item context effects.

# simulate data group 2
dat2 <- sirt::sim.raschtype( stats::rnorm(n2), b=b2, fixed.a=a2 )
colnames(dat2) <- paste0("I", 1:I, "des2" )
# define joint dataset
dat <- matrix( NA, nrow=n1+n2, ncol=2*I)
colnames(dat) <- c( colnames(dat1), colnames(dat2) )
dat[ 1:n1, 1:I ] <- dat1
dat[ n1 + 1:n2, I + 1:I ] <- dat2
# define group identifier
group <- c( rep(1,n1), rep(2,n2) )

#***
# Model 1: Rasch model two groups
itemindex <- rep( 1:I, 2 )
mod1 <- sirt::rasch.mml2( dat, group=group, est.b=itemindex )
summary(mod1)

#***
# Model 2: two item slope groups and designmatrix for intercepts
designmatrix <- matrix( 0, 2*I, I+1)
designmatrix[ ( 1:I )+ I,1:I] <- designmatrix[1:I,1:I] <- diag(I)
designmatrix[ ( 1:I )+ I,I+1] <- 1
mod2 <- sirt::rasch.mml2( dat, est.a=rep(1:2,each=I), designmatrix=designmatrix )
summary(mod2)

#############################################################################
# EXAMPLE 7: PIRLS dataset with missing responses
#############################################################################

data(data.pirlsmissing)
items <- grep( "R31", colnames(data.pirlsmissing), value=TRUE )
I <- length(items)
dat <- data.pirlsmissing

#****
# Model 1: recode missing responses as missing (missing are ignorable)

# data recoding
dat1 <- dat
dat1[ dat1==9 ] <- NA
# estimate Rasch model
mod1 <- sirt::rasch.mml2( dat1[,items], weights=dat$studwgt, group=dat$country )
summary(mod1)
##   Mean=0 0.341 -0.134 0.219
##   SD=1.142 1.166 1.197 0.959

#****
# Model 2: recode missing responses as wrong

# data recoding
dat2 <- dat
dat2[ dat2==9 ] <- 0
# estimate Rasch model
mod2 <- sirt::rasch.mml2( dat2[,items], weights=dat$studwgt, group=dat$country )
summary(mod2)
  ##   Mean=0 0.413 -0.172 0.446
  ##   SD=1.199 1.263 1.32 0.996

#****
# Model 3: recode missing responses as rho * P_i( theta ) and
#          apply pseudo-log-likelihood estimation
# Missing item responses are predicted by the model implied probability
# P_i( theta ) where theta is the ability estimate when ignoring missings (Model 1)
# and rho is an adjustment parameter. rho=0 is equivalent to Model 2 (treating
# missing as wrong) and rho=1 is equivalent to Model 1 (treating missing as ignorable).

# data recoding
dat3 <- dat
# simulate theta estimate from posterior distribution
theta <- stats::rnorm( nrow(dat3), mean=mod1$person$EAP, sd=mod1$person$SE.EAP )
rho <- .3   # define a rho parameter value of .3
for (ii in items){
    ind <- which( dat[,ii]==9 )
    dat3[ind,ii] <- rho*stats::plogis( theta[ind] - mod1$item$b[ which( items==ii ) ] )
                }

# estimate Rasch model
mod3 <- sirt::rasch.mml2( dat3[,items], weights=dat$studwgt, group=dat$country )
summary(mod3)
  ##   Mean=0 0.392 -0.153 0.38
  ##   SD=1.154 1.209 1.246 0.973

#****
# Model 4: simulate missing responses as rho * P_i( theta )
# The definition is the same as in Model 3. But it is now assumed
# that the missing responses are 'latent responses'.
set.seed(789)

# data recoding
dat4 <- dat
# simulate theta estimate from posterior distribution
theta <- stats::rnorm( nrow(dat4), mean=mod1$person$EAP, sd=mod1$person$SE.EAP )
rho <- .3   # define a rho parameter value of .3
for (ii in items){
    ind <- which( dat[,ii]==9 )
    p3 <- rho*stats::plogis( theta[ind] - mod1$item$b[ which( items==ii ) ] )
    dat4[ ind, ii ] <- 1*( stats::runif( length(ind), 0, 1 ) < p3)
                }

# estimate Rasch model
mod4 <- sirt::rasch.mml2( dat4[,items], weights=dat$studwgt, group=dat$country )
summary(mod4)
  ##   Mean=0 0.396 -0.156 0.382
  ##   SD=1.16 1.216 1.253 0.979

#****
# Model 5: recode missing responses for multiple choice items with four alternatives
#          to 1/4 and apply pseudo-log-likelihood estimation.
#          Missings for constructed response items are treated as incorrect.

# data recoding
dat5 <- dat
items_mc <- items[ substring( items, 7,7)=="M" ]
items_cr <- items[ substring( items, 7,7)=="C" ]
for (ii in items_mc){
    ind <- which( dat[,ii]==9 )
    dat5[ind,ii] <- 1/4
                }
for (ii in items_cr){
    ind <- which( dat[,ii]==9 )
    dat5[ind,ii] <- 0
                }

# estimate Rasch model
mod5 <- sirt::rasch.mml2( dat5[,items], weights=dat$studwgt, group=dat$country )
summary(mod5)
  ##   Mean=0 0.411 -0.165 0.435
  ##   SD=1.19 1.245 1.293 0.995

#*** For the following analyses, we ignore sample weights and the
#    country grouping.
data(data.pirlsmissing)
items <- grep( "R31", colnames(data.pirlsmissing), value=TRUE )
dat <- data.pirlsmissing
dat1 <- dat
dat1[ dat1==9 ] <- 0

#*** Model 6: estimate item difficulties assuming incorrect missing data treatment
mod6 <- sirt::rasch.mml2( dat1[,items], mmliter=50 )
summary(mod6)

#*** Model 7: reestimate model with constrained item difficulties
I <- length(items)
constraints <- cbind( 1:I, mod6$item$b )
mod7 <- sirt::rasch.mml2( dat1[,items], constraints=constraints)
summary(mod7)

#*** Model 8: score all missings responses as missing items
dat2 <- dat[,items]
dat2[ dat2==9 ] <- NA
mod8 <- sirt::rasch.mml2( dat2, constraints=constraints, mu.fixed=NULL )
summary(mod8)

#*** Model 9: estimate missing data model 'missing1' assuming a missingness
#       parameter delta.miss of zero
dat2 <-  dat[,items]    # note that missing item responses must be defined by 9
mod9 <- sirt::rasch.mml2( dat2, constraints=constraints, irtmodel="missing1",
            theta.k=seq(-5,5,len=10), delta.miss=0, mitermax=4, mu.fixed=NULL )
summary(mod9)

#*** Model 10: estimate missing data model with a large negative missing delta parameter
#=> This model is equivalent to treating missing responses as wrong
mod10 <- sirt::rasch.mml2( dat2, constraints=constraints, irtmodel="missing1",
             theta.k=seq(-5, 5, len=10), delta.miss=-10, mitermax=4, mmliter=200,
             mu.fixed=NULL )
summary(mod10)

#*** Model 11: choose a missingness delta parameter of -1
mod11 <- sirt::rasch.mml2( dat2, constraints=constraints, irtmodel="missing1",
             theta.k=seq(-5, 5, len=10), delta.miss=-1, mitermax=4,
             mmliter=200, mu.fixed=NULL )
summary(mod11)

#*** Model 12: estimate joint delta parameter
mod12 <- sirt::rasch.mml2( dat2, irtmodel="missing1", mu.fixed=cbind( c(1,2), 0 ),
             theta.k=seq(-8, 8, len=10), delta.miss=0, mitermax=4,
             mmliter=30, est.delta=rep(1,I)  )
summary(mod12)

#*** Model 13: estimate delta parameter in item groups defined by item format
est.delta <- 1 + 1 * ( substring( colnames(dat2),7,7 )=="M" )
mod13 <- sirt::rasch.mml2( dat2, irtmodel="missing1", mu.fixed=cbind( c(1,2), 0 ),
             theta.k=seq(-8, 8, len=10), delta.miss=0, mitermax=4,
             mmliter=30, est.delta=est.delta  )
summary(mod13)

#*** Model 14: estimate item specific delta parameter
mod14 <- sirt::rasch.mml2( dat2, irtmodel="missing1", mu.fixed=cbind( c(1,2), 0 ),
             theta.k=seq(-8, 8, len=10), delta.miss=0, mitermax=4,
             mmliter=30, est.delta=1:I  )
summary(mod14)

#############################################################################
# EXAMPLE 8: Comparison of different models for polytomous data
#############################################################################

data(data.Students, package="CDM")
head(data.Students)
dat <- data.Students[, paste0("act",1:5) ]
I <- ncol(dat)

#**************************************************
#*** Model 1: Partial Credit Model (PCM)

#*** Model 1a: PCM in TAM
mod1a <- TAM::tam.mml( dat )
summary(mod1a)

#*** Model 1b: PCM in sirt
mod1b <- sirt::rm.facets( dat )
summary(mod1b)

#*** Model 1c: PCM in mirt
mod1c <- mirt::mirt( dat, 1, itemtype=rep("Rasch",I), verbose=TRUE )
print(mod1c)

#**************************************************
#*** Model 2: Sequential Model (SM): Equal Loadings

#*** Model 2a: SM in sirt
dat1 <- CDM::sequential.items(dat)
resp <- dat1$dat.expand
iteminfo <- dat1$iteminfo
# fit model
mod2a <- sirt::rasch.mml2( resp )
summary(mod2a)

#**************************************************
#*** Model 3: Sequential Model (SM): Different Loadings

#*** Model 3a: SM in sirt
mod3a <- sirt::rasch.mml2( resp, est.a=iteminfo$itemindex )
summary(mod3a)

#**************************************************
#*** Model 4: Generalized partial credit model (GPCM)

#*** Model 4a: GPCM in TAM
mod4a <- TAM::tam.mml.2pl( dat, irtmodel="GPCM")
summary(mod4a)

#**************************************************
#*** Model 5: Graded response model (GRM)

#*** Model 5a: GRM in mirt
mod5a <- mirt::mirt( dat, 1, itemtype=rep("graded",I), verbose=TRUE)
print(mod5a)

# model comparison
logLik(mod1a);logLik(mod1b);mod1c@logLik  # PCM
logLik(mod2a)   # SM (Rasch)
logLik(mod3a)   # SM (GPCM)
logLik(mod4a)   # GPCM
mod5a@logLik    # GRM
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Marginal maximum likelihood (MML)}
%% \keyword{Ramsay's quotient model}
%% \keyword{Nonparametric item response theory}
%% \keyword{Pseudo-likelihood estimation}
%% \keyword{summary}
%% \keyword{plot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
