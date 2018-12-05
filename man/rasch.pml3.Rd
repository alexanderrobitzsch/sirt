%% File Name: rasch.pml3.Rd
%% File Version: 2.74



\name{rasch.pml3}
\alias{rasch.pml3}
\alias{summary.rasch.pml}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pairwise Marginal Likelihood Estimation for the Probit Rasch Model
}
\description{
This function estimates unidimensional 1PL and 2PL models with
the probit link using pairwise marginal maximum likelihood
estimation (PMML; Renard, Molenberghs & Geys, 2004).
Item pairs within an itemcluster can be excluded from the
pairwise likelihood (argument \code{itemcluster}).
The other alternative is to model a residual
error structure with itemclusters (argument \code{error.corr}).
}
\usage{
rasch.pml3(dat, est.b=seq(1, ncol(dat)), est.a=rep(0,ncol(dat)),
    est.sigma=TRUE, itemcluster=NULL, weight=rep(1, nrow(dat)), numdiff.parm=0.001,
    b.init=NULL, a.init=NULL,  sigma.init=NULL, error.corr=0*diag( 1, ncol(dat) ),
    err.constraintM=NULL, err.constraintV=NULL, glob.conv=10^(-6), conv1=10^(-4),
    pmliter=300, progress=TRUE, use.maxincrement=TRUE )

\method{summary}{rasch.pml}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
An \eqn{N \times I} data frame of dichotomous item responses
}
  \item{est.b}{
Vector of integers of length \eqn{I}. Same integers mean that the
corresponding items do have the same item difficulty \code{b}.
Entries of \code{0} mean fixing item parameters to values
specified in \code{b.init}.
}
\item{est.a}{
Vector of integers of length \eqn{I}. Same integers mean that the
corresponding items do have the same item slope \code{a}.
Entries of \code{0} mean fixing item parameters to values
specified in \code{a.init}.
}
  \item{est.sigma}{
Should sigma (the trait standard deviation) be estimated?
The default is \code{TRUE}.
}
  \item{itemcluster}{
Optional vector of length \eqn{I} of integers which indicates itemclusters.
Same integers correspond to the same itemcluster. An entry of \code{0}
correspond to an item which is not included in any itemcluster.
}
  \item{weight}{
Optional vector of person weights
}
  \item{numdiff.parm}{
Step parameter for numerical differentiation
}
  \item{b.init}{
Initial or fixed item difficulty
}
  \item{a.init}{
Initial or fixed item slopes
}
  \item{sigma.init}{
Initial or fixed trait standard deviation
}
  \item{error.corr}{
An optional \eqn{I \times I} integer matrix
which defines the estimation of residual correlations.
Entries of zero indicate that the corresponding
residual correlation should not be estimated.
Integers which differ from zero indicate correlations to be estimated.
All entries with an equal integer are estimated by the same residual
correlation. The default of \code{error.corr} is a diagonal matrix
which means that no residual correlation is estimated. If \code{error.corr}
deviates from this default, then the argument \code{itemcluster} is set
to \code{NULL}.\cr
If some error correlations are estimated, then no itempairs in
\code{itemcluster} can be excluded from the pairwise modeling.
}
%     err.constraintM=NULL, err.constraintV=NULL,
\item{err.constraintM}{
An optional \eqn{P \times L} matrix where \eqn{P} denotes
the number of item pairs in pseudolikelihood estimation
and \eqn{L} is the number of linear constraints for residual
correlations (see Details).
        }
\item{err.constraintV}{
An optional \eqn{L \times 1} matrix with specified values
for linear constraints on residual correlations (see Details).
        }
  \item{glob.conv}{
Global convergence criterion
}
  \item{conv1}{
Convergence criterion for model parameters
}
  \item{pmliter}{
Maximum number of iterations
}
  \item{progress}{
Display progress?
}
\item{use.maxincrement}{Optional logical whether increments in
    slope parameters should be controlled in size in iterations.
    The default is \code{TRUE}.
}
  \item{object}{
Object of class \code{rasch.pml}
}
  \item{\dots}{
Further arguments to be passed
}
}
\details{
The probit item response model can be estimated with this function:
\deqn{P(X_{pi}=1|\theta_p )=\Phi( a_i  \theta_p - b_i ) \quad, \quad
    \theta_p \sim N ( 0, \sigma^2 )}{%
P(X_pi=1|\theta_p )=\Phi( a_i  \theta_p - b_i ), \theta_p ~ N ( 0, \sigma^2 )
    }
where \eqn{\Phi} denotes the normal distribution function. This model can
also be expressed as a latent variable model which assumes
a latent response tendency \eqn{X_{pi}^\ast}{X_pi*} which is equal to
1 if \eqn{X_{pi}> - b_i} and otherwise zero. If \eqn{\epsilon_{pi}} is
standard normally distributed, then
\deqn{X_{pi}^{\ast}=a_i \theta_p - b_i + \epsilon_{pi} }{%
 X_pi*=a_i \theta_p - b_i + \epsilon_{pi}
            }
An arbitrary pattern of residual correlations between
\eqn{\epsilon_{pi}} and \eqn{\epsilon_{pj}} for item pairs \eqn{i}
and \eqn{j} can be imposed using the \code{error.corr} argument.

Linear constraints \eqn{Me=v} on residual correlations
\eqn{e=Cov( \epsilon_{pi}, \epsilon_{pj})_{ij}} (in a vectorized form) can be
specified using the arguments \code{err.constraintM} (matrix \eqn{M})
and \code{err.constraintV} (vector \eqn{v}). The estimation
is described in Neuhaus (1996).

For the pseudo likelihood information criterion (PLIC)
see Stanford and Raftery (2002).
}
\value{
A list with following entries:
\item{item}{Data frame with estimated item parameters}
\item{iter}{Number of iterations}
\item{deviance}{Pseudolikelihood multiplied by minus 2}
\item{b}{Estimated item difficulties}
\item{sigma}{Estimated standard deviation}
\item{dat}{Original dataset}
\item{ic}{Data frame with information criteria (sample size,
number of estimated parameters, pseudolikelihood
information criterion \code{PLIC})}
\item{link}{Used link function (only probit is permitted)}
\item{itempairs}{Estimated statistics of item pairs}
\item{error.corr}{Estimated error correlation matrix}
\item{eps.corr}{
    Vectorized error correlation matrix
        }
\item{omega.rel}{
    Reliability of the sum score according to Green and Yang (2009).
    If some item pairs are excluded in the estimation, the residual
    correlation for these item pairs is assumed to be zero.
        }

\item{\dots}{}
}

\references{
Green, S. B., & Yang, Y. (2009). Reliability of summed item
scores using structural equation modeling: An alternative to
coefficient alpha. \emph{Psychometrika, 74}, 155-167.

Neuhaus, W. (1996). Optimal estimation under
linear constraints. \emph{Astin Bulletin, 26}, 233-245.

Renard, D., Molenberghs, G., & Geys, H. (2004). A pairwise likelihood
approach to estimation in multilevel probit models. \emph{Computational Statistics
& Data Analysis, 44}, 649-667.

Stanford, D. C., & Raftery, A. E. (2002).
Approximate Bayes factors for image segmentation:
The pseudolikelihood information criterion (PLIC).
\emph{IEEE Transactions on Pattern Analysis and
Machine Intelligence, 24}, 1517-1520.
}




%\author{
%Alexander Robitzsch
%}


\note{
This function needs the \pkg{combinat} library.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Get a summary of \code{rasch.pml3} with \code{\link{summary.rasch.pml}}.

For simulation of locally dependent items see \code{\link{sim.rasch.dep}}.

For pairwise conditional likelihood estimation see \code{\link{rasch.pairwise}}
or \code{\link{rasch.pairwise.itemcluster}}.

For an assessment of global model fit see \code{\link{modelfit.sirt}}.
}

%% print % operator in example by using \%
%% see the tam.mml manual
%% a <- E \%*\% lambda


\examples{
#############################################################################
# EXAMPLE 1: Reading data set
#############################################################################

data(data.read)
dat <- data.read

#******
# Model 1: Rasch model with PML estimation
mod1 <- sirt::rasch.pml3( dat )
summary(mod1)

#******
# Model 2: Excluding item pairs with local dependence
#          from bivariate composite likelihood
itemcluster <- rep( 1:3, each=4)
mod2 <- sirt::rasch.pml3( dat, itemcluster=itemcluster )
summary(mod2)

\dontrun{
#*****
# Model 3: Modelling error correlations:
#          joint residual correlations for each itemcluster
error.corr <- diag(1,ncol(dat))
for ( ii in 1:3){
    ind.ii <- which( itemcluster==ii )
    error.corr[ ind.ii, ind.ii ] <- ii
        }
# estimate the model with error correlations
mod3 <- sirt::rasch.pml3( dat, error.corr=error.corr )
summary(mod3)

#****
# Model 4: model separate residual correlations
I <- ncol(error.corr)
error.corr1 <- matrix( 1:(I*I), ncol=I )
error.corr <- error.corr1 * ( error.corr > 0 )
# estimate the model with error correlations
mod4 <- sirt::rasch.pml3( dat, error.corr=error.corr )
summary(mod4)

#****
# Model 5:  assume equal item difficulties:
# b_1=b_7 and b_2=b_12
# fix item difficulty of the 6th item to .1
est.b <- 1:I
est.b[7] <- 1; est.b[12] <- 2 ; est.b[6] <- 0
b.init <- rep( 0, I ) ; b.init[6] <- .1
mod5 <- sirt::rasch.pml3( dat, est.b=est.b, b.init=b.init)
summary(mod5)

#****
# Model 6: estimate three item slope groups
est.a <- rep(1:3, each=4 )
mod6 <- sirt::rasch.pml3( dat, est.a=est.a, est.sigma=0)
summary(mod6)

#############################################################################
# EXAMPLE 2: PISA reading
#############################################################################

data(data.pisaRead)
dat <- data.pisaRead$data

# select items
dat <- dat[, substring(colnames(dat),1,1)=="R" ]

#******
# Model 1: Rasch model with PML estimation
mod1 <- sirt::rasch.pml3( as.matrix(dat) )
  ## Trait SD (Logit Link) : 1.419

#******
# Model 2: Model correlations within testlets
error.corr <- diag(1,ncol(dat))
testlets <- paste( data.pisaRead$item$testlet )
itemcluster <- match( testlets, unique(testlets ) )
for ( ii in 1:(length(unique(testlets))) ){
    ind.ii <- which( itemcluster==ii )
    error.corr[ ind.ii, ind.ii ] <- ii
        }
# estimate the model with error correlations
mod2 <- sirt::rasch.pml3( dat, error.corr=error.corr )
  ## Trait SD (Logit Link) : 1.384

#****
# Model 3: model separate residual correlations
I <- ncol(error.corr)
error.corr1 <- matrix( 1:(I*I), ncol=I )
error.corr <- error.corr1 * ( error.corr > 0 )
# estimate the model with error correlations
mod3 <- sirt::rasch.pml3( dat, error.corr=error.corr )
  ## Trait SD (Logit Link) : 1.384

#############################################################################
# EXAMPLE 3: 10 locally independent items
#############################################################################

#**********
# simulate some data
set.seed(554)
N <- 500    # persons
I <- 10        # items
theta <- stats::rnorm(N,sd=1.3 )    # trait SD of 1.3
b <- seq(-2, 2, length=I) # item difficulties

# simulate data from the Rasch model
dat <- sirt::sim.raschtype( theta=theta, b=b )

# estimation with rasch.pml and probit link
mod1 <- sirt::rasch.pml3( dat )
summary(mod1)

# estimation with rasch.mml2 function
mod2 <- sirt::rasch.mml2( dat )

# estimate item parameters for groups with five item parameters each
est.b <- rep( 1:(I/2), each=2 )
mod3 <- sirt::rasch.pml3( dat, est.b=est.b )
summary(mod3)

# compare parameter estimates
summary(mod1)
summary(mod2)
summary(mod3)

#############################################################################
# EXAMPLE 4: 11 items and 2 item clusters with 2 and 3 items
#############################################################################

set.seed(5698)
I <- 11                             # number of items
n <- 5000                           # number of persons
b <- seq(-2,2, len=I)               # item difficulties
theta <- stats::rnorm( n, sd=1 ) # person abilities
# itemcluster
itemcluster <- rep(0,I)
itemcluster[c(3,5)] <- 1
itemcluster[c(2,4,9)] <- 2
# residual correlations
rho <- c( .7, .5 )

# simulate data (under the logit link)
dat <- sirt::sim.rasch.dep( theta, b, itemcluster, rho )
colnames(dat) <- paste("I", seq(1,ncol(dat)), sep="")

#***
# Model 1: estimation using the Rasch model (with probit link)
mod1 <- sirt::rasch.pml3( dat )
#***
# Model 2: estimation when pairs of locally dependent items are eliminated
mod2 <- sirt::rasch.pml3( dat, itemcluster=itemcluster)

#***
# Model 3: Positive correlations within testlets
est.corrs <- diag( 1, I )
est.corrs[ c(3,5), c(3,5) ] <- 2
est.corrs[ c(2,4,9), c(2,4,9) ] <- 3
mod3 <- sirt::rasch.pml3( dat, error.corr=est.corrs )

#***
# Model 4: Negative correlations between testlets
est.corrs <- diag( 1, I )
est.corrs[ c(3,5), c(2,4,9) ] <- 2
est.corrs[ c(2,4,9), c(3,5) ] <- 2
mod4 <- sirt::rasch.pml3( dat, error.corr=est.corrs )

#***
# Model 5: sum constraint of zero within and between testlets
est.corrs <- matrix( 1:(I*I),  I, I )
cluster2 <- c(2,4,9)
est.corrs[ setdiff( 1:I, c(cluster2)),  ] <- 0
est.corrs[, setdiff( 1:I, c(cluster2))  ] <- 0
# define an error constraint matrix
itempairs0 <- mod4$itempairs
IP <- nrow(itempairs0)
err.constraint <- matrix( 0, IP, 1 )
err.constraint[ ( itempairs0$item1 \%in\% cluster2 )
       & ( itempairs0$item2 \%in\% cluster2 ), 1 ] <- 1
# set sum of error covariances to 1.2
err.constraintV <- matrix(3*.4,1,1)

mod5 <- sirt::rasch.pml3( dat, error.corr=est.corrs,
         err.constraintM=err.constraint, err.constraintV=err.constraintV)

#****
# Model 6: Constraint on sum of all correlations
est.corrs <- matrix( 1:(I*I),  I, I )
# define an error constraint matrix
itempairs0 <- mod4$itempairs
IP <- nrow(itempairs0)
# define two side conditions
err.constraint <- matrix( 0, IP, 2 )
err.constraintV <- matrix( 0, 2, 1)
# sum of all correlations is zero
err.constraint[, 1 ] <- 1
err.constraintV[1,1] <- 0
# sum of items cluster c(1,2,3) is 0
cluster2 <- c(1,2,3)
err.constraint[ ( itempairs0$item1 \%in\%  cluster2 )
       & ( itempairs0$item2 \%in\% cluster2 ), 2 ] <- 1
err.constraintV[2,1] <- 0

mod6 <- sirt::rasch.pml3( dat, error.corr=est.corrs,
    err.constraintM=err.constraint,  err.constraintV=err.constraintV)
summary(mod6)

#############################################################################
# EXAMPLE 5: 10 Items: Cluster 1 -> Items 1,2
#         Cluster 2 -> Items 3,4,5;   Cluster 3 -> Items 7,8,9
#############################################################################

set.seed(7650)
I <- 10                             # number of items
n <- 5000                           # number of persons
b <- seq(-2,2, len=I)               # item difficulties
bsamp <- b <- sample(b)             # sample item difficulties
theta <- stats::rnorm( n, sd=1 ) # person abilities
# define itemcluster
itemcluster <- rep(0,I)
itemcluster[ 1:2 ] <- 1
itemcluster[ 3:5 ] <- 2
itemcluster[ 7:9 ] <- 3
# define residual correlations
rho <- c( .55, .35, .45)

# simulate data
dat <- sirt::sim.rasch.dep( theta, b, itemcluster, rho )
colnames(dat) <- paste("I", seq(1,ncol(dat)), sep="")

#***
# Model 1: residual correlation (equal within item clusters)
# define a matrix of integers for estimating error correlations
error.corr <- diag(1,ncol(dat))
for ( ii in 1:3){
    ind.ii <- which( itemcluster==ii )
    error.corr[ ind.ii, ind.ii ] <- ii
        }
# estimate the model
mod1 <- sirt::rasch.pml3( dat, error.corr=error.corr )

#***
# Model 2: residual correlation (different within item clusters)
# define again a matrix of integers for estimating error correlations
error.corr <- diag(1,ncol(dat))
for ( ii in 1:3){
    ind.ii <- which( itemcluster==ii )
    error.corr[ ind.ii, ind.ii ] <- ii
        }
I <- ncol(error.corr)
error.corr1 <- matrix( 1:(I*I), ncol=I )
error.corr <- error.corr1 * ( error.corr > 0 )
# estimate the model
mod2 <- sirt::rasch.pml3( dat, error.corr=error.corr )

#***
# Model 3: eliminate item pairs within itemclusters for PML estimation
mod3 <- sirt::rasch.pml3( dat, itemcluster=itemcluster )

#***
# Model 4: Rasch model ignoring dependency
mod4 <- sirt::rasch.pml3( dat )

# compare different models
summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Pairwise marginal maximum likelihood (PMML)}
\keyword{Local dependence}
\keyword{summary}
