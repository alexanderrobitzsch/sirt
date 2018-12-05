%% File Name: rasch.mirtlc.Rd
%% File Version: 2.932

\name{rasch.mirtlc}
\alias{rasch.mirtlc}
\alias{summary.rasch.mirtlc}
\alias{logLik.rasch.mirtlc}
\alias{anova.rasch.mirtlc}
\alias{IRT.irfprob.rasch.mirtlc}
\alias{IRT.likelihood.rasch.mirtlc}
\alias{IRT.posterior.rasch.mirtlc}
\alias{IRT.modelfit.rasch.mirtlc}
\alias{summary.IRT.modelfit.rasch.mirtlc}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multidimensional Latent Class 1PL and 2PL Model
}
\description{
This function estimates the multidimensional latent class
Rasch (1PL) and 2PL model (Bartolucci, 2007; Bartolucci, Montanari & Pandolfi,
2012) for dichotomous data which emerges from the original latent class model
(Goodman, 1974) and a multidimensional IRT model.
}
\usage{
rasch.mirtlc(dat, Nclasses=NULL, modeltype="LC", dimensions=NULL,
    group=NULL, weights=rep(1,nrow(dat)), theta.k=NULL, ref.item=NULL,
    distribution.trait=FALSE,  range.b=c(-8,8), range.a=c(.2, 6 ),
    progress=TRUE, glob.conv=10^(-5), conv1=10^(-5), mmliter=1000,
    mstep.maxit=3, seed=0, nstarts=1, fac.iter=.35)

\method{summary}{rasch.mirtlc}(object,...)

\method{anova}{rasch.mirtlc}(object,...)

\method{logLik}{rasch.mirtlc}(object,...)

\method{IRT.irfprob}{rasch.mirtlc}(object,...)

\method{IRT.likelihood}{rasch.mirtlc}(object,...)

\method{IRT.posterior}{rasch.mirtlc}(object,...)

\method{IRT.modelfit}{rasch.mirtlc}(object,...)

\method{summary}{IRT.modelfit.rasch.mirtlc}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
An \eqn{N \times I} data frame
}
  \item{Nclasses}{
Number of latent classes. If the trait vector (or matrix)
\code{theta.k} is specified, then \code{Nclasses} is set
to the dimension of \code{theta.k}.
}
  \item{modeltype}{
Modeltype. \code{LC} is the latent class model of Goodman (1974).
\code{MLC1} is the multidimensional latent class Rasch model with
item discrimination parameter of 1. \code{MLC2} allows for the
estimation of item discriminations.
}
\item{dimensions}{
    Vector of dimension integers which allocate items to dimensions.
        }
  \item{group}{
A group identifier for multiple group estimation
}
  \item{weights}{
 Vector of sample weights
 }
  \item{theta.k}{
A grid of theta values can be specified if theta
should not be estimated. In the one-dimensional case, it
must be a vector, in the \eqn{D}-dimensional case it must
be a matrix of dimension \eqn{D}.
}
\item{ref.item}{An optional vector of integers which indicate the items
    whose intercept and slope are fixed at 0 and 1, respectively.
    }
\item{distribution.trait}{A type of the assumed theta distribution can
be specified. One alternative is \code{normal}
for the normal distribution assumption.
The options \code{smooth2}, \code{smooth3} and
\code{smooth4} use the log-linear smoothing of
 Xu and von Davier (2008) to smooth the distribution
 up to two, three or four moments, respectively.
This function only works in unidimensional models. \cr
If a different string is provided as an input (e.g. \code{no}),
then no smoothing is conducted.
    }
  \item{range.b}{
 Range of item difficulties which are allowed for estimation
 }
  \item{range.a}{
 Range of item slopes which are allowed for estimation
 }
  \item{progress}{
Display progress? Default is \code{TRUE}.
}
  \item{glob.conv}{
Global relative deviance convergence criterion
}
  \item{conv1}{
Item parameter convergence criterion
}
  \item{mmliter}{
Maximum number of iterations
}
  \item{mstep.maxit}{
Maximum number of iterations within an M step
}
  \item{seed}{
Set random seed for latent class estimation. A seed
can be specified. If the seed is negative, then
the function will generate a random seed.
}
\item{nstarts}{
    If a positive integer is provided, then a \code{nstarts}
    starts with different starting values are conducted.
}
\item{fac.iter}{A parameter between 0 and 1 to control the maximum
increment in each iteration. The larger the parameter the more increments
will become smaller from iteration to iteration.
}
  \item{object}{
Object of class \code{rasch.mirtlc}
}
  \item{\dots}{
Further arguments to be passed
}
}
\details{
The multidimensional latent class Rasch model (Bartolucci, 2007)
is an item response model which combines ideas from
latent class analysis and item response models with continuous variables.
With \code{modeltype="MLC2"} the following \eqn{D}-dimensional
item response model is estimated
\deqn{logit P(X_{pi}=1 | \theta_p )=a_i \theta_{pcd}- b_i}{%
logit P(X_{pi}=1 | \theta_p )=a_i \theta_{pcd}- b_i }
Besides the item thresholds \eqn{b_i} and item slopes \eqn{a_i},
for a prespecified number of latent classes \eqn{c=1,\ldots,C}{c=1,...,C}
a set of \eqn{C} \eqn{D}-dimensional \eqn{\{\theta_{cd} \}_{cd}}{{\theta_{cd} }_{cd}}
vectors are estimated.
These vectors represent the locations of latent classes. If the user
provides a grid of theta distribution \code{theta.k} as an argument in
\code{rasch.mirtlc}, then the ability distribution is fixed.

In the unidimensional Rasch model with \eqn{I} items, \eqn{(I+1)/2}
(if \eqn{I} odd) or \eqn{I/2 + 1} (if \eqn{I} even) trait location
parameters are identified (see De Leeuw & Verhelst, 1986;
Lindsay et al., 1991; for a review see Formann, 2007).
}


\value{
A list with following entries
\item{pjk}{Item probabilities evaluated at discretized ability distribution}
\item{rprobs}{Item response probabilities like in \code{pjk},
    but for each item category}
\item{pi.k}{Estimated trait distribution}
\item{theta.k}{Discretized ability distribution}
\item{item}{Estimated item parameters}
\item{trait}{Estimated ability distribution (\code{theta.k} and \code{pi.k})}
\item{mean.trait}{Estimated mean of ability distribution}
\item{sd.trait}{Estimated standard deviation of ability distribution}
\item{skewness.trait}{Estimated skewness of ability distribution}
\item{cor.trait}{Estimated correlation between abilities (only applies for
    multidimensional models)}
\item{ic}{Information criteria}
\item{D}{Number of dimensions}
\item{G}{Number of groups}
\item{deviance}{Deviance}
\item{ll}{Log-likelihood}
\item{Nclasses}{Number of classes}
\item{modeltype}{Used model type}
\item{estep.res}{Result from E step: \code{f.qk.yi} is the individual posterior,
    \code{f.yi.qk} is the individual likelihood}
\item{dat}{Original data frame}
\item{devL}{Vector of deviances if multiple random starts were conducted}
\item{seedL}{Vector of seed if multiple random starts were conducted}
\item{iter}{Number of iterations}
}



\references{
Bartolucci, F. (2007). A class of multidimensional
IRT models for testing unidimensionality and clustering
items. \emph{Psychometrika, 72}, 141-157.

Bartolucci, F., Montanari, G. E., & Pandolfi, S. (2012).
Dimensionality of the latent structure and item selection via latent
class multidimensional IRT models.
\emph{Psychometrika, 77}, 782-802.

De Leeuw, J., & Verhelst, N. (1986). Maximum likelihood estimation in
generalized Rasch models. \emph{Journal of Educational and Behavioral Statistics,
11}, 183-196.

Formann, A. K. (2007). (Almost) Equivalence between conditional and mixture
maximum likelihood estimates for some models of the Rasch type.
In M. von Davier & C. H. Carstensen:
\emph{Multivariate and Mixture Distribution Rasch Models} (pp. 177-189).
Springer: New York.

Goodman, L. A. (1974). Exploratory latent structure
analysis using both identifiable and unidentifiable
models. \emph{Biometrika, 61}, 215-231.

Lindsay, B., Clogg, C. C., & Grego, J. (1991).
Semiparametric estimation in the Rasch model and related exponential response
models, including a simple latent class model for item analysis.
\emph{Journal of the American Statistical Association, 86}, 96-107.

Xu, X., & von Davier, M. (2008). \emph{Fitting the structured general
diagnostic model to NAEP data}. ETS Research Report ETS RR-08-27.
Princeton, ETS.
}


%\author{
%Alexander Robitzsch
%}



\note{
For the estimation of latent class models,
rerunning the model with different starting values
(different random seeds) is recommended.

For fixed theta estimation in the multidimensional
case, large vectors are generated during estimation
leading to memory overflow in \R.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \link[CDM]{\code{gdm}}
\seealso{
See also the \code{\link[CDM:gdm]{CDM::gdm}}
function in the \pkg{CDM} package.

For an assessment of global model fit see \code{\link{modelfit.sirt}}.

The estimation of the multidimensional latent class item response
model for polytomous data can be conducted in the \pkg{MultiLCIRT}
package. Latent class analysis can be carried out with
\pkg{poLCA} and \pkg{randomLCA} packages.
}
\examples{
#############################################################################
# EXAMPLE 1: Reading data
#############################################################################
data( data.read )
dat <- data.read

#***************
# latent class models

# latent class model with 1 class
mod1 <- sirt::rasch.mirtlc( dat, Nclasses=1 )
summary(mod1)

# latent class model with 2 classes
mod2 <- sirt::rasch.mirtlc( dat, Nclasses=2 )
summary(mod2)

\dontrun{
# latent class model with 3 classes
mod3 <- sirt::rasch.mirtlc( dat, Nclasses=3, seed=- 30)
summary(mod3)

# extract individual likelihood
lmod3 <- IRT.likelihood(mod3)
str(lmod3)
# extract likelihood value
logLik(mod3)
# extract item response functions
IRT.irfprob(mod3)

# compare models 1, 2 and 3
anova(mod2,mod3)
IRT.compareModels(mod1,mod2,mod3)
# avsolute and relative model fit
smod2 <- IRT.modelfit(mod2)
smod3 <- IRT.modelfit(mod3)
summary(smod2)
IRT.compareModels(smod2,smod3)

# latent class model with 4 classes and 3 starts with different seeds
mod4 <- sirt::rasch.mirtlc( dat, Nclasses=4,seed=-30,  nstarts=3 )
# display different solutions
sort(mod4$devL)
summary(mod4)

# latent class multiple group model
# define group identifier
group <- rep( 1, nrow(dat))
group[ 1:150 ] <- 2
mod5 <- sirt::rasch.mirtlc( dat, Nclasses=3, group=group )
summary(mod5)

#*************
# Unidimensional IRT models with ordered trait

# 1PL model with 3 classes
mod11 <- sirt::rasch.mirtlc( dat, Nclasses=3, modeltype="MLC1", mmliter=30)
summary(mod11)

# 1PL model with 11 classes
mod12 <- sirt::rasch.mirtlc( dat, Nclasses=11,modeltype="MLC1", mmliter=30)
summary(mod12)

# 1PL model with 11 classes and fixed specified theta values
mod13 <- sirt::rasch.mirtlc( dat,  modeltype="MLC1",
             theta.k=seq( -4, 4, len=11 ), mmliter=100)
summary(mod13)

# 1PL model with fixed theta values and normal distribution
mod14 <- sirt::rasch.mirtlc( dat,  modeltype="MLC1", mmliter=30,
             theta.k=seq( -4, 4, len=11 ), distribution.trait="normal")
summary(mod14)

# 1PL model with a smoothed trait distribution (up to 3 moments)
mod15 <- sirt::rasch.mirtlc( dat,  modeltype="MLC1", mmliter=30,
             theta.k=seq( -4, 4, len=11 ),  distribution.trait="smooth3")
summary(mod15)

# 2PL with 3 classes
mod16 <- sirt::rasch.mirtlc( dat, Nclasses=3, modeltype="MLC2", mmliter=30 )
summary(mod16)

# 2PL with fixed theta and smoothed distribution
mod17 <- sirt::rasch.mirtlc( dat, theta.k=seq(-4,4,len=12), mmliter=30,
             modeltype="MLC2", distribution.trait="smooth4"  )
summary(mod17)

# 1PL multiple group model with 8 classes
# define group identifier
group <- rep( 1, nrow(dat))
group[ 1:150 ] <- 2
mod21 <- sirt::rasch.mirtlc( dat, Nclasses=8, modeltype="MLC1", group=group )
summary(mod21)

#***************
# multidimensional latent class IRT models

# define vector of dimensions
dimensions <- rep( 1:3, each=4 )

# 3-dimensional model with 8 classes and seed 145
mod31 <- sirt::rasch.mirtlc( dat, Nclasses=8, mmliter=30,
             modeltype="MLC1", seed=145, dimensions=dimensions )
summary(mod31)

# try the model above with different starting values
mod31s <- sirt::rasch.mirtlc( dat, Nclasses=8,
             modeltype="MLC1", seed=-30, nstarts=30, dimensions=dimensions )
summary(mod31s)

# estimation with fixed theta vectors
#=> 4^3=216 classes
theta.k <- seq(-4, 4, len=6 )
theta.k <- as.matrix( expand.grid( theta.k, theta.k, theta.k ) )
mod32 <- sirt::rasch.mirtlc( dat,  dimensions=dimensions,
              theta.k=theta.k, modeltype="MLC1"  )
summary(mod32)

# 3-dimensional 2PL model
mod33 <- sirt::rasch.mirtlc( dat, dimensions=dimensions, theta.k=theta.k, modeltype="MLC2")
summary(mod33)

#############################################################################
# EXAMPLE 2: Skew trait distribution
#############################################################################
set.seed(789)
N <- 1000   # number of persons
I <- 20     # number of items
theta <- sqrt( exp( stats::rnorm( N ) ) )
theta <- theta - mean(theta )
# calculate skewness of theta distribution
mean( theta^3 ) / stats::sd(theta)^3
# simulate item responses
dat <- sirt::sim.raschtype( theta, b=seq(-2,2,len=I ) )

# normal distribution
mod1 <- sirt::rasch.mirtlc( dat, theta.k=seq(-4,4,len=15), modeltype="MLC1",
               distribution.trait="normal", mmliter=30)

# allow for skew distribution with smoothed distribution
mod2 <- sirt::rasch.mirtlc( dat, theta.k=seq(-4,4,len=15), modeltype="MLC1",
               distribution.trait="smooth3", mmliter=30)

# nonparametric distribution
mod3 <- sirt::rasch.mirtlc( dat, theta.k=seq(-4,4,len=15), modeltype="MLC1", mmliter=30)

summary(mod1)
summary(mod2)
summary(mod3)

#############################################################################
# EXAMPLE 3: Stouffer-Toby dataset data.si02 with 5 items
#############################################################################

data(dat.si02)
dat <- data.si02$data
weights <- data.si02$weights   # extract weights

# Model 1: 2 classes Rasch model
mod1 <- sirt::rasch.mirtlc( dat, Nclasses=2, modeltype="MLC1", weights=weights,
                 ref.item=4, nstarts=5)
summary(mod1)

# Model 2: 3 classes Rasch model: not all parameters are identified
mod2 <- sirt::rasch.mirtlc( dat, Nclasses=3, modeltype="MLC1", weights=weights,
                ref.item=4, nstarts=5)
summary(mod2)

# Model 3: Latent class model with 2 classes
mod3 <- sirt::rasch.mirtlc( dat, Nclasses=2, modeltype="LC", weights=weights, nstarts=5)
summary(mod3)

# Model 4: Rasch model with normal distribution
mod4 <- sirt::rasch.mirtlc( dat,  modeltype="MLC1", weights=weights,
            theta.k=seq( -6, 6, len=21 ), distribution.trait="normal", ref.item=4)
summary(mod4)
}

#############################################################################
# EXAMPLE 4: 5 classes, 3 dimensions and 27 items
#############################################################################

set.seed(979)
I <- 9
N <- 5000
b <- seq( - 1.5, 1.5, len=I)
b <- rep(b,3)
# define class locations
theta.k <- c(-3.0, -4.1, -2.8, 1.7, 2.3, 1.8,
   0.2, 0.4, -0.1,   2.6, 0.1, -0.9, -1.1,-0.7, 0.9 )

Nclasses <- 5
theta.k0 <- theta.k <- matrix( theta.k, Nclasses, 3, byrow=TRUE )
pi.k <- c(.20,.25,.25,.10,.15)
theta <- theta.k[ rep( 1:Nclasses, round(N*pi.k) ), ]
dimensions <- rep( 1:3, each=I)
# simulate item responses
dat <- matrix( NA, nrow=N, ncol=I*3)
for (ii in 1:(3*I) ){
    dat[,ii] <- 1 * ( stats::runif(N) < stats::plogis( theta[,dimensions[ii]] - b[ii]))
}
colnames(dat) <- paste0( rep( LETTERS[1:3], each=I ), 1:(3*I) )

# estimate model
mod1 <- sirt::rasch.mirtlc( dat, Nclasses=Nclasses, dimensions=dimensions,
             modeltype="MLC1", ref.item=c(5,14,23), glob.conv=.0005, conv1=.0005)

round( cbind( mod1$theta.k, mod1$pi.k ), 3 )
  ##          [,1]   [,2]   [,3]  [,4]
  ##   [1,] -2.776 -3.791 -2.667 0.250
  ##   [2,] -0.989 -0.605  0.957 0.151
  ##   [3,]  0.332  0.418 -0.046 0.246
  ##   [4,]  2.601  0.171 -0.854 0.101
  ##   [5,]  1.791  2.330  1.844 0.252
cbind( theta.k, pi.k )
  ##                       pi.k
  ##   [1,] -3.0 -4.1 -2.8 0.20
  ##   [2,]  1.7  2.3  1.8 0.25
  ##   [3,]  0.2  0.4 -0.1 0.25
  ##   [4,]  2.6  0.1 -0.9 0.10
  ##   [5,] -1.1 -0.7  0.9 0.15

# plot class locations
plot( 1:3, mod1$theta.k[1,], xlim=c(1,3), ylim=c(-5,3), col=1, pch=1, type="n",
    axes=FALSE, xlab="Dimension", ylab="Location")
axis(1, 1:3 ) ;  axis(2) ; axis(4)
for (cc in 1:Nclasses){ # cc <- 1
    lines(1:3, mod1$theta.k[cc,], col=cc, lty=cc )
    points(1:3, mod1$theta.k[cc,], col=cc,  pch=cc )
}

\dontrun{
#------
# estimate model with gdm function in CDM package
library(CDM)
# define Q-matrix
Qmatrix <- matrix(0,3*I,3)
Qmatrix[ cbind( 1:(3*I), rep(1:3, each=I) ) ] <- 1

set.seed(9176)
# random starting values for theta locations
theta.k <- matrix( 2*stats::rnorm(5*3), 5, 3 )
colnames(theta.k) <- c("Dim1","Dim2","Dim3")
# try possibly different starting values

# estimate model in CDM
b.constraint  <- cbind( c(5,14,23), 1, 0 )
mod2 <- CDM::gdm( dat, theta.k=theta.k, b.constraint=b.constraint, skillspace="est",
               irtmodel="1PL",  Qmatrix=Qmatrix)
summary(mod2)

#------
# estimate model with MultiLCIRT package
miceadds::library_install("MultiLCIRT")

# define matrix to allocate each item to one dimension
multi1 <- matrix( 1:(3*I), nrow=3, byrow=TRUE )
# define reference items in item-dimension allocation matrix
multi1[ 1, c(1,5)  ] <- c(5,1)
multi1[ 2, c(10,14) - 9  ] <- c(14,9)
multi1[ 3, c(19,23) - 18 ] <- c(23,19)

# Rasch model with 5 latent classes (random start: start=1)
mod3 <- MultiLCIRT::est_multi_poly(S=dat,k=5,       # k=5 ability levels
                start=1,link=1,multi=multi1,tol=10^-5,
                output=TRUE, disp=TRUE, fort=TRUE)
# estimated location points and class probabilities in MultiLCIRT
cbind( t( mod3$Th ), mod3$piv )
# compare results with rasch.mirtlc
cbind( mod1$theta.k, mod1$pi.k )
# simulated data parameters
cbind( theta.k, pi.k )

#----
# estimate model with cutomized input in mirt
library(mirt)
#-- define Theta design matrix for 5 classes
Theta <- diag(5)
Theta <- cbind( Theta, Theta, Theta )
r1 <- rownames(Theta) <- paste0("C",1:5)
colnames(Theta) <- c( paste0(r1, "D1"), paste0(r1, "D2"), paste0(r1, "D3") )
  ##      C1D1 C2D1 C3D1 C4D1 C5D1 C1D2 C2D2 C3D2 C4D2 C5D2 C1D3 C2D3 C3D3 C4D3 C5D3
  ##   C1    1    0    0    0    0    1    0    0    0    0    1    0    0    0    0
  ##   C2    0    1    0    0    0    0    1    0    0    0    0    1    0    0    0
  ##   C3    0    0    1    0    0    0    0    1    0    0    0    0    1    0    0
  ##   C4    0    0    0    1    0    0    0    0    1    0    0    0    0    1    0
  ##   C5    0    0    0    0    1    0    0    0    0    1    0    0    0    0    1
#-- define mirt model
I <- ncol(dat)  # I=27
mirtmodel <- mirt::mirt.model("
        C1D1=1-9 \n C2D1=1-9 \n  C3D1=1-9 \n  C4D1=1-9  \n  C5D1=1-9
        C1D2=10-18 \n C2D2=10-18 \n  C3D2=10-18 \n  C4D2=10-18  \n  C5D2=10-18
        C1D3=19-27 \n C2D3=19-27 \n  C3D3=19-27 \n  C4D3=19-27  \n  C5D3=19-27
        CONSTRAIN=(1-9,a1),(1-9,a2),(1-9,a3),(1-9,a4),(1-9,a5),
                    (10-18,a6),(10-18,a7),(10-18,a8),(10-18,a9),(10-18,a10),
                    (19-27,a11),(19-27,a12),(19-27,a13),(19-27,a14),(19-27,a15)
                ")
#-- get initial parameter values
mod.pars <- mirt::mirt(dat, model=mirtmodel,  pars="values")
#-- redefine initial parameter values
# set all d parameters initially to zero
ind <- which( ( mod.pars$name=="d" ) )
mod.pars[ ind,"value" ]  <- 0
# fix item difficulties of reference items to zero
mod.pars[ ind[ c(5,14,23) ], "est"] <- FALSE
mod.pars[ind,]
# initial item parameters of cluster locations (a1,...,a15)
ind <- which( ( mod.pars$name \%in\% paste0("a", c(1,6,11) ) ) & ( mod.pars$est ) )
mod.pars[ind,"value"] <- -2
ind <- which( ( mod.pars$name \%in\% paste0("a", c(1,6,11)+1 ) ) & ( mod.pars$est ) )
mod.pars[ind,"value"] <- -1
ind <- which( ( mod.pars$name \%in\% paste0("a", c(1,6,11)+2 ) ) & ( mod.pars$est ) )
mod.pars[ind,"value"] <- 0
ind <- which( ( mod.pars$name \%in\% paste0("a", c(1,6,11)+3 ) ) & ( mod.pars$est ) )
mod.pars[ind,"value"] <- 1
ind <- which( ( mod.pars$name \%in\% paste0("a", c(1,6,11)+4 ) ) & ( mod.pars$est ) )
mod.pars[ind,"value"] <- 0
#-- define prior for latent class analysis
lca_prior <- function(Theta,Etable){
  TP <- nrow(Theta)
  if ( is.null(Etable) ){ prior <- rep( 1/TP, TP ) }
  if ( ! is.null(Etable) ){
    prior <- ( rowSums(Etable[, seq(1,2*I,2)]) + rowSums(Etable[,seq(2,2*I,2)]) )/I
  }
  prior <- prior / sum(prior)
  return(prior)
}

#-- estimate model in mirt
mod4 <- mirt::mirt(dat, mirtmodel, pars=mod.pars, verbose=TRUE,
              technical=list( customTheta=Theta, customPriorFun=lca_prior,
                    MAXQUAD=1E20) )
# correct number of estimated parameters
mod4@nest <- as.integer(sum(mod.pars$est) + nrow(Theta)-1 )
# extract coefficients
# source.all(pfsirt)
cmod4 <- sirt::mirt.wrapper.coef(mod4)

# estimated item difficulties
dfr <- data.frame( "sim"=b, "mirt"=-cmod4$coef$d, "sirt"=mod1$item$thresh )
round( dfr, 4 )
  ##         sim    mirt    sirt
  ##   1  -1.500 -1.3782 -1.3382
  ##   2  -1.125 -1.0059 -0.9774
  ##   3  -0.750 -0.6157 -0.6016
  ##   4  -0.375 -0.2099 -0.2060
  ##   5   0.000  0.0000  0.0000
  ##   6   0.375  0.5085  0.4984
  ##   7   0.750  0.8661  0.8504
  ##   8   1.125  1.3079  1.2847
  ##   9   1.500  1.5891  1.5620
  ##   [...]

#-- reordering estimated latent clusters to make solutions comparable
#* extract estimated cluster locations from sirt
order.sirt <- c(1,5,3,4,2)  # sort(order.sirt)
round(mod1$trait[,1:3],3)
dfr <- data.frame( "sim"=theta.k, mod1$trait[order.sirt,1:3] )
colnames(dfr)[4:6] <- paste0("sirt",1:3)
#* extract estimated cluster locations from mirt
c4 <- cmod4$coef[, paste0("a",1:15) ]
c4 <- apply( c4,2, FUN=function(ll){ ll[ ll!=0 ][1] } )
trait.loc <- matrix(c4,5,3)
order.mirt <- c(1,4,3,5,2)  # sort(order.mirt)
dfr <- cbind( dfr, trait.loc[ order.mirt, ] )
colnames(dfr)[7:9] <- paste0("mirt",1:3)
# compare estimated cluster locations
round(dfr,3)
  ##     sim.1 sim.2 sim.3  sirt1  sirt2  sirt3  mirt1  mirt2  mirt3
  ##   1  -3.0  -4.1  -2.8 -2.776 -3.791 -2.667 -2.856 -4.023 -2.741
  ##   5   1.7   2.3   1.8  1.791  2.330  1.844  1.817  2.373  1.869
  ##   3   0.2   0.4  -0.1  0.332  0.418 -0.046  0.349  0.421 -0.051
  ##   4   2.6   0.1  -0.9  2.601  0.171 -0.854  2.695  0.166 -0.876
  ##   2  -1.1  -0.7   0.9 -0.989 -0.605  0.957 -1.009 -0.618  0.962
#* compare estimated cluster sizes
dfr <- data.frame( "sim"=pi.k, "sirt"=mod1$pi.k[order.sirt,1],
            "mirt"=mod4@Prior[[1]][ order.mirt] )
round(dfr,4)
  ##      sim   sirt   mirt
  ##   1 0.20 0.2502 0.2500
  ##   2 0.25 0.2522 0.2511
  ##   3 0.25 0.2458 0.2494
  ##   4 0.10 0.1011 0.0986
  ##   5 0.15 0.1507 0.1509

#############################################################################
# EXAMPLE 5: Dataset data.si04 from Bartolucci et al. (2012)
#############################################################################

data(data.si04)

# define reference items
ref.item <- c(7,17,25,44,64)
dimensions <- data.si04$itempars$dim

# estimate a Rasch latent class with 9 classes
mod1 <- sirt::rasch.mirtlc( data.si04$data, Nclasses=9, dimensions=dimensions,
             modeltype="MLC1", ref.item=ref.item, glob.conv=.005, conv1=.005,
             nstarts=1, mmliter=200 )

# compare estimated distribution with simulated distribution
round( cbind( mod1$theta.k, mod1$pi.k ), 4 ) # estimated
  ##            [,1]    [,2]    [,3]    [,4]    [,5]   [,6]
  ##    [1,] -3.6043 -5.1323 -5.3022 -6.8255 -4.3611 0.1341
  ##    [2,]  0.2083 -2.7422 -2.8754 -5.3416 -2.5085 0.1573
  ##    [3,] -2.8641 -4.0272 -5.0580 -0.0340 -0.9113 0.1163
  ##    [4,] -0.3575 -2.0081 -1.7431  1.2992 -0.1616 0.0751
  ##    [5,]  2.9329  0.3662 -1.6516 -3.0284  0.1844 0.1285
  ##    [6,]  1.5092 -2.0461 -4.3093  1.0481  1.0806 0.1094
  ##    [7,]  3.9899  3.1955 -4.0010  1.8879  2.2988 0.1460
  ##    [8,]  4.3062  0.7080 -1.2324  1.4351  2.0893 0.1332
  ##    [9,]  5.0855  4.1214 -0.9141  2.2744  1.5314 0.0000

round(d2,4) # simulated
  ##         class      A      B      C      D      E     pi
  ##    [1,]     1 -3.832 -5.399 -5.793 -7.042 -4.511 0.1323
  ##    [2,]     2 -2.899 -4.217 -5.310 -0.055 -0.915 0.1162
  ##    [3,]     3 -0.376 -2.137 -1.847  1.273 -0.078 0.0752
  ##    [4,]     4  0.208 -2.934 -3.011 -5.526 -2.511 0.1583
  ##    [5,]     5  1.536 -2.137 -4.606  1.045  1.143 0.1092
  ##    [6,]     6  2.042 -0.573 -0.404 -4.331 -1.044 0.0471
  ##    [7,]     7  3.853  0.841 -2.993 -2.746  0.803 0.0822
  ##    [8,]     8  4.204  3.296 -4.328  1.892  2.419 0.1453
  ##    [9,]     9  4.466  0.700 -1.334  1.439  2.161 0.1343
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Multidimensional latent class Rasch model}
\keyword{summary}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
