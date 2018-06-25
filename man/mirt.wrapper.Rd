%% File Name: mirt.wrapper.Rd
%% File Version: 1.81

\name{mirt.wrapper}
\alias{mirt.wrapper}
\alias{mirt.wrapper.coef}
\alias{mirt.wrapper.posterior}
\alias{mirt.wrapper.fscores}
\alias{mirt.wrapper.itemplot}
\alias{IRT.irfprob.SingleGroupClass}
\alias{IRT.likelihood.SingleGroupClass}
\alias{IRT.posterior.SingleGroupClass}

% SingleGroupClass

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Some Functions for Wrapping with the \pkg{mirt} Package
}
\description{
Some functions for wrapping with the \pkg{mirt} package.
}
\usage{
# extract coefficients
mirt.wrapper.coef(mirt.obj)

# extract posterior, likelihood, ...
mirt.wrapper.posterior(mirt.obj,weights=NULL)
\method{IRT.likelihood}{SingleGroupClass}(object, \dots)
\method{IRT.posterior}{SingleGroupClass}(object, \dots)

# S3 method for extracting item response functions
\method{IRT.irfprob}{SingleGroupClass}(object, \dots)

# compute factor scores
mirt.wrapper.fscores(mirt.obj,weights=NULL)

# convenience function for itemplot
mirt.wrapper.itemplot( mirt.obj, ask=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{mirt.obj}{
A fitted model in \pkg{mirt} package
}
\item{object}{A fitted object in \pkg{mirt} package of class
    \code{SingleGroupClass} or \code{MultipleGroupClass}.
}
\item{weights}{Optional vector of student weights}
\item{ask}{Optional logical indicating whether each new plot should be
confirmed.}
\item{\dots}{Further arguments to be passed.}
}

\details{
The function \code{mirt.wrapper.coef} collects all item parameters
in a data frame.

The function \code{mirt.wrapper.posterior} extracts the individual
likelihood, individual likelihood and expected counts. This function does not
yet cover the case of multiple groups.

The function \code{mirt.wrapper.fscores} computes factor scores
EAP, MAP and MLE. The factor scores are computed on the
discrete grid of latent traits (contrary to the computation in \code{mirt}) as
specified in \code{mirt.obj@Theta}. This function does also not work
for multiple groups.

The function \code{mirt.wrapper.itemplot} displays all item plots
after each other.
}

\value{

Function \code{mirt.wrapper.coef} -- List with entries

\item{coef}{Data frame with item parameters}
\item{GroupPars}{Data frame or list with distribution parameters}


Function \code{mirt.wrapper.posterior} -- List with entries

\item{theta.k}{Grid of theta points}
\item{pi.k}{Trait distribution on \code{theta.k}}
\item{f.yi.qk}{Individual likelihood}
\item{f.qk.yi}{Individual posterior}
\item{n.ik}{Expected counts}
\item{data}{Used dataset}

Function \code{mirt.wrapper.fscores} -- List with entries

\item{person}{Data frame with person parameter estimates (factor scores)
    EAP, MAP and MLE for all dimensions.
}
\item{EAP.rel}{EAP reliabilities}
}


%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Alexander Robitzsch
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\seealso{
See the \pkg{mirt} package on CRAN
\url{https://CRAN.R-project.org/package=mirt}
and on GitHub \url{https://github.com/philchalmers/mirt}.
% \href{https://github.com/philchalmers/mirt}{GitHub}

See \url{https://groups.google.com/forum/#!forum/mirt-package}
for discussion about the \pkg{mirt} package. \cr


See for the main estimation functions in \pkg{mirt}:
\code{\link[mirt:mirt]{mirt::mirt}},
\code{\link[mirt:multipleGroup]{mirt::multipleGroup}}
and \code{\link[mirt:bfactor]{mirt::bfactor}}.

See \code{\link[mirt:coef-method]{mirt::coef-method}}for extracting
coefficients.

See \code{\link[mirt:mod2values]{mirt::mod2values}} for collecting
parameter values in a mirt parameter table.

See \code{\link{lavaan2mirt}} for converting \code{lavaan} syntax
to \code{mirt} syntax.

See \code{\link{tam2mirt}} for converting fitted \code{tam} models
into \code{mirt} objects.

See also \code{\link[CDM:IRT.likelihood]{CDM::IRT.likelihood}},
\code{\link[CDM:IRT.posterior]{CDM::IRT.posterior}} and
\code{\link[CDM:IRT.irfprob]{CDM::IRT.irfprob}} for general
extractor functions.
}


% See \code{\link{mirt.model.vars}} for specifying \pkg{mirt} models
% using variable names.

%% ~~objects to See Also as \code{\link{help}}, ~~~




\section{Examples for the \pkg{mirt} Package}{

%\itemize{
\enumerate{

%1
\item Latent class analysis (\code{\link{data.read}}, Model 7)

%2
\item Mixed Rasch model  (\code{\link{data.read}}, Model 8)

%3
\item Located unidimensional and multidimensional
latent class models / Multidimensional latent class IRT models
(\code{\link{data.read}}, Model 12;
\code{\link{rasch.mirtlc}}, Example 4)

%4
\item Multidimensional IRT model with discrete latent traits
 (\code{\link{data.read}}, Model 13)

%5
\item DINA model (\code{\link{data.read}}, Model 14;
\code{\link[CDM:data.dcm]{data.dcm}}, \pkg{CDM}, Model 1m)

%6
\item Unidimensional IRT model with non-normal distribution
 (\code{\link{data.read}}, Model 15)

%7
\item Grade of membership model  (\code{\link{gom.em}}, Example 2)

%8
\item Rasch copula model (\code{\link{rasch.copula2}}, Example 5)

%9
\item Additive GDINA model
(\code{\link[CDM:data.dcm]{data.dcm}}, \pkg{CDM}, Model 6m)

%10
\item Longitudinal Rasch model (\code{\link{data.long}}, Model 3)

%11
\item Normally distributed residuals (\code{\link{data.big5}}, Example 1, Model 5)

% 12
\item Nedelsky model (\code{\link{nedelsky.irf}},
Examples 1, 2)

% 13
\item Beta item response model (\code{\link{brm.irf}}, Example 1)
}

}


\examples{
\dontrun{
# A development version can be installed from GitHub
if (FALSE){ # default is set to FALSE, use the installed version
   library(devtools)
   devtools::install_github("philchalmers/mirt")
          }
# now, load mirt
library(mirt)

#############################################################################
# EXAMPLE 1: Extracting item parameters and posterior LSAT data
#############################################################################

data(LSAT7, package="mirt")
data <- mirt::expand.table(LSAT7)

#*** Model 1: 3PL model for item 5 only, other items 2PL
mod1 <- mirt::mirt(data, 1, itemtype=c("2PL","2PL","2PL","2PL","3PL"), verbose=TRUE)
print(mod1)
summary(mod1)
# extracting coefficients
coef(mod1)
mirt.wrapper.coef(mod1)$coef
# extract parameter values in mirt
mirt::mod2values(mod1)
# extract posterior
post1 <- sirt::mirt.wrapper.posterior(mod1)
# extract item response functions
probs1 <- IRT.irfprob( mod1 )
str(probs1)
# extract individual likelihood
likemod1 <- IRT.likelihood( mod1 )
str(likemod1)
# extract individual posterior
postmod1 <- IRT.posterior( mod1 )
str(postmod1)

#*** Model 2: Confirmatory model with two factors
cmodel <- mirt::mirt.model("
        F1=1,4,5
        F2=2,3
        ")
mod2 <- mirt::mirt(data, cmodel, verbose=TRUE)
print(mod2)
summary(mod2)
# extract coefficients
coef(mod2)
mirt.wrapper.coef(mod2)$coef
# extract posterior
post2 <- sirt::mirt.wrapper.posterior(mod2)

#############################################################################
# EXAMPLE 2: Extracting item parameters and posterior for differering
#            number of response catagories | Dataset Science
#############################################################################

data(Science,package="mirt")
library(psych)
psych::describe(Science)

# modify dataset
dat <- Science
dat[ dat[,1] > 3,1] <- 3
psych::describe(dat)

# estimate generalized partial credit model
mod1 <- mirt::mirt(dat, 1, itemtype="gpcm")
print(mod1)
# extract coefficients
coef(mod1)
mirt.wrapper.coef(mod1)$coef
# extract posterior
post1 <- sirt::mirt.wrapper.posterior(mod1)

#############################################################################
# EXAMPLE 3: Multiple group model; simulated dataset from mirt package
#############################################################################

#*** simulate data (copy from the multipleGroup manual site in mirt package)
set.seed(1234)
a <- matrix(c(abs( stats::rnorm(5,1,.3)), rep(0,15),abs( stats::rnorm(5,1,.3)),
          rep(0,15),abs( stats::rnorm(5,1,.3))), 15, 3)
d <- matrix( stats::rnorm(15,0,.7),ncol=1)
mu <- c(-.4, -.7, .1)
sigma <- matrix(c(1.21,.297,1.232,.297,.81,.252,1.232,.252,1.96),3,3)
itemtype <- rep("dich", nrow(a))
N <- 1000
dataset1 <- mirt::simdata(a, d, N, itemtype)
dataset2 <- mirt::simdata(a, d, N, itemtype, mu=mu, sigma=sigma)
dat <- rbind(dataset1, dataset2)
group <- c(rep("D1", N), rep("D2", N))

#group models
model <- mirt::mirt.model("
   F1=1-5
   F2=6-10
   F3=11-15
      ")

# separate analysis
mod_configural <- mirt::multipleGroup(dat, model, group=group, verbose=TRUE)
mirt.wrapper.coef(mod_configural)

# equal slopes (metric invariance)
mod_metric <- mirt::multipleGroup(dat, model, group=group, invariance=c("slopes"),
                verbose=TRUE)
mirt.wrapper.coef(mod_metric)

# equal slopes and intercepts (scalar invariance)
mod_scalar <- mirt::multipleGroup(dat, model, group=group,
          invariance=c("slopes","intercepts","free_means","free_varcov"), verbose=TRUE)
mirt.wrapper.coef(mod_scalar)

# full constraint
mod_fullconstrain <- mirt::multipleGroup(dat, model, group=group,
             invariance=c("slopes", "intercepts"), verbose=TRUE )
mirt.wrapper.coef(mod_fullconstrain)

#############################################################################
# EXAMPLE 4: Nonlinear item response model
#############################################################################

data(data.read)
dat <- data.read
# specify mirt model with some interactions
mirtmodel <- mirt.model("
   A=1-4
   B=5-8
   C=9-12
   (A*B)=4,8
   (C*C)=9
   (A*B*C)=12
   " )
# estimate model
res <- mirt::mirt( dat, mirtmodel, verbose=TRUE, technical=list(NCYCLES=3) )
# look at estimated parameters
mirt.wrapper.coef(res)
coef(res)
mirt::mod2values(res)
# model specification
res@model

#############################################################################
# EXAMPLE 5: Extracting factor scores
#############################################################################

data(data.read)
dat <- data.read
# define lavaan model and convert syntax to mirt
lavmodel <- "
    A=~ a*A1+a*A2+1.3*A3+A4       # set loading of A3 to 1.3
    B=~ B1+1*B2+b3*B3+B4
    C=~ c*C1+C2+c*C3+C4
    A1 | da*t1
    A3 | da*t1
    C4 | dg*t1
    B1 | 0*t1
    B3 | -1.4*t1                  # fix item threshold of B3 to -1.4
    A ~~ B                        # estimate covariance between A and B
    A ~~ .6 * C                   # fix covariance to .6
    B ~~ B                        # estimate variance of B
    A ~ .5*1                      # set mean of A to .5
    B ~ 1                         # estimate mean of B
    "
res <- sirt::lavaan2mirt( dat, lavmodel, verbose=TRUE, technical=list(NCYCLES=3) )
# estimated coefficients
mirt.wrapper.coef(res$mirt)
# extract factor scores
fres <- sirt::mirt.wrapper.fscores(res$mirt)
# look at factor scores
head( round(fres$person,2))
  ##     case    M EAP.Var1 SE.EAP.Var1 EAP.Var2 SE.EAP.Var2 EAP.Var3 SE.EAP.Var3 MLE.Var1
  ##   1    1 0.92     1.26        0.67     1.61        0.60     0.05        0.69     2.65
  ##   2    2 0.58     0.06        0.59     1.14        0.55    -0.80        0.56     0.00
  ##   3    3 0.83     0.86        0.66     1.15        0.55     0.48        0.74     0.53
  ##   4    4 1.00     1.52        0.67     1.57        0.60     0.73        0.76     2.65
  ##   5    5 0.50    -0.13        0.58     0.85        0.48    -0.82        0.55    -0.53
  ##   6    6 0.75     0.41        0.63     1.09        0.54     0.27        0.71     0.00
  ##     MLE.Var2 MLE.Var3 MAP.Var1 MAP.Var2 MAP.Var3
  ##   1     2.65    -0.53     1.06     1.59     0.00
  ##   2     1.06    -1.06     0.00     1.06    -1.06
  ##   3     1.06     2.65     1.06     1.06     0.53
  ##   4     2.65     2.65     1.59     1.59     0.53
  ##   5     0.53    -1.06    -0.53     0.53    -1.06
  ##   6     1.06     2.65     0.53     1.06     0.00
# EAP reliabilities
round(fres$EAP.rel,3)
  ##    Var1  Var2  Var3
  ##   0.574 0.452 0.541
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{mirt}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
