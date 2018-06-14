%% File Name: rasch.copula.Rd
%% File Version: 1.62

\name{rasch.copula2}
\alias{rasch.copula2}
\alias{rasch.copula3}
\alias{summary.rasch.copula2}
\alias{summary.rasch.copula3}
\alias{anova.rasch.copula2}
\alias{anova.rasch.copula3}
\alias{logLik.rasch.copula2}
\alias{logLik.rasch.copula3}
\alias{IRT.likelihood.rasch.copula2}
\alias{IRT.likelihood.rasch.copula3}
\alias{IRT.posterior.rasch.copula2}
\alias{IRT.posterior.rasch.copula3}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multidimensional IRT Copula Model
}
\description{
This function handles local dependence by specifying
copulas for residuals in multidimensional item response
models for dichotomous item responses
(Braeken, 2011; Braeken, Tuerlinckx & de Boeck, 2007;
Schroeders, Robitzsch & Schipolowski, 2014).
Estimation is allowed for item difficulties, item slopes and
a generalized logistic link function (Stukel, 1988).

The function \code{rasch.copula3} allows the estimation of multidimensional
models while \code{rasch.copula2} only handles unidimensional models.
}
\usage{
rasch.copula2(dat, itemcluster, copula.type ="bound.mixt" ,
    progress = TRUE, mmliter = 1000, delta = NULL,
    theta.k = seq(-4, 4, len = 21), alpha1 = 0, alpha2 = 0,
    numdiff.parm = 1e-06,  est.b = seq(1, ncol(dat)),
    est.a = rep(1, ncol(dat)), est.delta = NULL, b.init = NULL , a.init = NULL ,
    est.alpha = FALSE, glob.conv = 0.0001, alpha.conv = 1e-04, conv1 = 0.001,
    dev.crit=.2 , increment.factor=1.01)

rasch.copula3(dat, itemcluster, dims=NULL , copula.type ="bound.mixt" ,
    progress = TRUE, mmliter = 1000, delta = NULL,
    theta.k = seq(-4, 4, len = 21), alpha1 = 0, alpha2 = 0,
    numdiff.parm = 1e-06,  est.b = seq(1, ncol(dat)),
    est.a = rep(1, ncol(dat)), est.delta = NULL, b.init = NULL , a.init = NULL ,
    est.alpha = FALSE, glob.conv = 0.0001, alpha.conv = 1e-04, conv1 = 0.001,
    dev.crit=.2 , rho.init=.5 , increment.factor=1.01)

\method{summary}{rasch.copula2}(object,...)
\method{summary}{rasch.copula3}(object,...)

\method{anova}{rasch.copula2}(object,...)
\method{anova}{rasch.copula3}(object,...)

\method{logLik}{rasch.copula2}(object,...)
\method{logLik}{rasch.copula3}(object,...)

\method{IRT.likelihood}{rasch.copula2}(object,...)
\method{IRT.likelihood}{rasch.copula3}(object,...)

\method{IRT.posterior}{rasch.copula2}(object,...)
\method{IRT.posterior}{rasch.copula3}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
An \eqn{N \times I} data frame. Cases with only missing responses
are removed from the analysis.
}
  \item{itemcluster}{
An integer vector of length \eqn{I} (number of items). Items with
the same integers define a joint item cluster of (positively) locally
dependent items. Values of zero indicate that the corresponding item
is not included in any item cluster of dependent responses.
}
\item{dims}{
A vector indicating to which dimension an item is allocated.
The default is that all items load on the first dimension.
}
\item{copula.type}{
A character or a vector containing one of the following copula
types: \code{bound.mixt} (boundary mixture copula),
\code{cook.johnson} (Cook-Johnson copula) or \code{frank} (Frank copula)
(see Braeken, 2011).
The vector \code{copula.type} must match the number of different
itemclusters. For every itemcluster, a different copula type
may be specified (see Examples).
}  %% elaborate this further in the Details section
  \item{progress}{
Print progress? Default is \code{TRUE}.
}
  \item{mmliter}{
Maximum number of iterations.
}
  \item{delta}{
An optional vector of starting values for the dependency parameter \code{delta}.
}
  \item{theta.k}{
Discretized trait distribution
}
  \item{alpha1}{
\code{alpha1} parameter in the generalized logistic item response model
(Stukel, 1988). The default is 0 which leads together with \code{alpha2=0}
to the logistic link function.
}
  \item{alpha2}{
\code{alpha2} parameter in the generalized logistic item response model
}
  \item{numdiff.parm}{
Parameter for numerical differentiation
}
  \item{est.b}{
Integer vector of item difficulties to be estimated
}
  \item{est.a}{
Integer vector of item discriminations to be estimated
}
  \item{est.delta}{
Integer vector of length \code{length(itemcluster)}. Nonzero integers
correspond to \code{delta} parameters which are estimated.
Equal integers indicate parameter equality constraints.
}
  \item{b.init}{
Initial \eqn{b} parameters
}
  \item{a.init}{
Initial \eqn{a} parameters
}
  \item{est.alpha}{
Should both alpha parameters be estimated? Default is \code{FALSE}.
}
 \item{glob.conv}{
Convergence criterion for all parameters
}
  \item{alpha.conv}{
Maximal change in alpha parameters for convergence
}
  \item{conv1}{
Maximal change in item parameters for convergence
}
  \item{dev.crit}{
Maximal change in the deviance. Default is \code{.2}.
}
\item{rho.init}{Initial value for off-diagonal elements in correlation matrix}
\item{increment.factor}{A numeric value larger than one which controls the
    size of increments in iterations. To stabilize convergence,
    choose values 1.05 or 1.1 in some situations.}
  \item{object}{
Object of class \code{rasch.copula2} or \code{rasch.copula3}
}
  \item{\dots}{
Further arguments to be passed
}
}
%\details{
%% Details about alpha1 and alpha2 for generalized logistic function?
%}
\value{
A list with following entries
\item{N.itemclusters}{Number of item clusters}
\item{item}{Estimated item parameters}
\item{iter}{Number of iterations}
\item{dev}{Deviance}
\item{delta}{Estimated dependency parameters \eqn{\delta}}
\item{b}{Estimated item difficulties}
\item{a}{Estimated item slopes}
\item{mu}{Mean}
\item{sigma}{Standard deviation}
\item{alpha1}{Parameter \eqn{\alpha_1} in the generalized item response model}
\item{alpha2}{Parameter \eqn{\alpha_2} in the generalized item response model}
\item{ic}{Information criteria}
\item{theta.k}{Discretized ability distribution}
\item{pi.k}{Fixed \eqn{\theta} distribution}
\item{deviance}{Deviance}
\item{pattern}{Item response patterns with frequencies and posterior
distribution}
\item{person}{Data frame with person parameters}
\item{datalist}{List of generated data frames during estimation}
\item{EAP.rel}{Reliability of the EAP}
\item{copula.type}{Type of copula}
\item{summary.delta}{Summary for estimated \eqn{\delta} parameters}
\item{f.qk.yi}{Individual posterior}
\item{f.yi.qk}{Individual likelihood}
\item{\dots}{}
}
\references{
    Braeken, J. (2011). A boundary mixture approach to violations
    of conditional independence. \emph{Psychometrika}, \bold{76}, 57-76.

Braeken, J., & Tuerlinckx, F. (2009). Investigating latent constructs with item
response models: A MATLAB IRTm toolbox.
\emph{Behavior Research Methods}, \bold{41}, 1127-1137.

    Braeken, J., Tuerlinckx, F., & De Boeck, P. (2007).
    Copulas for residual dependencies.
    \emph{Psychometrika}, \bold{72}, 393-411.

Schroeders, U., Robitzsch, A., & Schipolowski, S. (2014). A comparison of different
psychometric approaches to modeling testlet structures: An example with C-tests.
\emph{Journal of Educational Measurement, 51}(4), 400-418.

    Stukel, T. A. (1988). Generalized logistic models.
    \emph{Journal of the American Statistical Association},
    \bold{83}, 426-431.
}
\author{
Alexander Robitzsch
}
%\note{
%The function \code{rasch.copula2} is substantially faster than
%\code{rasch.copula}, but it is in a more experimental status.
%}

%\section{Warning}{
%    The multiple group option only works for \code{rasch.copula},
%    but not for \code{rasch.copula2}. This option will be implemented
%    in the near future.
%                }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For a summary see \code{\link{summary.rasch.copula2}}.

For simulating locally dependent item responses see \code{\link{sim.rasch.dep}}.

Person parameters estimates are obtained by \code{\link{person.parameter.rasch.copula}}.

See \code{\link{rasch.mml2}} for the generalized logistic link function.

See also Braeken and Tuerlinckx (2009) for alternative (and more expanded)
copula models implemented in the MATLAB software. See
\url{http://ppw.kuleuven.be/okp/software/irtm/}.
}


\examples{
#############################################################################
# EXAMPLE 1: Reading Data
#############################################################################

data(data.read)
dat <- data.read

# define item clusters
itemcluster <- rep( 1:3 , each=4 )

# estimate Copula model
mod1 <- sirt::rasch.copula2( dat=dat , itemcluster=itemcluster)

# estimate Rasch model
mod2 <- sirt::rasch.copula2( dat=dat , itemcluster=itemcluster ,
        delta=rep(0,3) , est.delta=rep(0,3 ) )
summary(mod1)
summary(mod2)

\dontrun{
#############################################################################
# EXAMPLE 2: 11 items nested within 2 item clusters (testlets)
#    with 2 resp. 3 dependent and 6 independent items
#############################################################################

set.seed(5698)
I <- 11                             # number of items
n <- 3000                           # number of persons
b <- seq(-2,2, len=I)               # item difficulties
theta <- stats::rnorm( n , sd = 1 ) # person abilities
# define item clusters
itemcluster <- rep(0,I)
itemcluster[ c(3,5 )] <- 1
itemcluster[c(2,4,9)] <- 2
# residual correlations
rho <- c( .7 , .5 )

# simulate data
dat <- sirt::sim.rasch.dep( theta , b , itemcluster , rho )
colnames(dat) <- paste("I" , seq(1,ncol(dat)) , sep="")

# estimate Rasch copula model
mod1 <- sirt::rasch.copula2( dat , itemcluster = itemcluster )
summary(mod1)

# both item clusters have Cook-Johnson copula as dependency
mod1c <- sirt::rasch.copula2( dat , itemcluster = itemcluster ,
            copula.type ="cook.johnson")
summary(mod1c)

# first item boundary mixture and second item Cook-Johnson copula
mod1d <- sirt::rasch.copula2( dat , itemcluster = itemcluster ,
            copula.type = c( "bound.mixt" , "cook.johnson" ) )
summary(mod1d)

# compare result with Rasch model estimation in rasch.copula2
# delta must be set to zero
mod2 <- sirt::rasch.copula2( dat , itemcluster = itemcluster , delta = c(0,0) ,
            est.delta = c(0,0) )
summary(mod2)

#############################################################################
# EXAMPLE 3: 12 items nested within 3 item clusters (testlets)
#   Cluster 1 -> Items 1-4; Cluster 2 -> Items 6-9;  Cluster 3 -> Items 10-12
#############################################################################

set.seed(967)
I <- 12                             # number of items
n <- 450                            # number of persons
b <- seq(-2,2, len=I)               # item difficulties
b <- sample(b)                      # sample item difficulties
theta <- stats::rnorm( n , sd = 1 ) # person abilities
# itemcluster
itemcluster <- rep(0,I)
itemcluster[ 1:4 ] <- 1
itemcluster[ 6:9 ] <- 2
itemcluster[ 10:12 ] <- 3
# residual correlations
rho <- c( .35 , .25 , .30 )

# simulate data
dat <- sirt::sim.rasch.dep( theta , b , itemcluster , rho )
colnames(dat) <- paste("I" , seq(1,ncol(dat)) , sep="")

# estimate Rasch copula model
mod1 <- sirt::rasch.copula2( dat , itemcluster = itemcluster )
summary(mod1)

# person parameter estimation assuming the Rasch copula model
pmod1 <- sirt::person.parameter.rasch.copula(raschcopula.object = mod1 )

# Rasch model estimation
mod2 <- sirt::rasch.copula2( dat , itemcluster = itemcluster ,
             delta = rep(0,3) , est.delta = rep(0,3) )
summary(mod1)
summary(mod2)

#############################################################################
# EXAMPLE 4: Two-dimensional copula model
#############################################################################

set.seed(5698)
I <- 9
n <- 1500                           # number of persons
b <- seq(-2,2, len=I)               # item difficulties
theta0 <- stats::rnorm( n , sd = sqrt( .6 ) )

#*** Dimension 1
theta <- theta0 + stats::rnorm( n , sd = sqrt( .4 ) )   # person abilities
# itemcluster
itemcluster <- rep(0,I)
itemcluster[ c(3,5 )] <- 1
itemcluster[c(2,4,9)] <- 2
itemcluster1 <- itemcluster
# residual correlations
rho <- c( .7 , .5 )
# simulate data
dat <- sirt::sim.rasch.dep( theta , b , itemcluster , rho )
colnames(dat) <- paste("A" , seq(1,ncol(dat)) , sep="")
dat1 <- dat
# estimate model of dimension 1
mod0a <- sirt::rasch.copula2( dat1 , itemcluster = itemcluster1)
summary(mod0a)

#*** Dimension 2
theta <- theta0 + stats::rnorm( n , sd = sqrt( .8 ) )        # person abilities
# itemcluster
itemcluster <- rep(0,I)
itemcluster[ c(3,7,8 )] <- 1
itemcluster[c(4,6)] <- 2
itemcluster2 <- itemcluster
# residual correlations
rho <- c( .2, .4 )
# simulate data
dat <- sirt::sim.rasch.dep( theta , b , itemcluster , rho )
colnames(dat) <- paste("B" , seq(1,ncol(dat)) , sep="")
dat2 <- dat
# estimate model of dimension 2
mod0b <- sirt::rasch.copula2( dat2 , itemcluster = itemcluster2)
summary(mod0b)

# both dimensions
dat <- cbind( dat1 , dat2 )
itemcluster2 <- ifelse( itemcluster2 > 0 , itemcluster2 +2 , 0 )
itemcluster <- c( itemcluster1 , itemcluster2 )
dims <- rep( 1:2 , each=I)

# estimate two-dimensional copula model
mod1 <- sirt::rasch.copula3( dat , itemcluster=itemcluster , dims=dims , est.a=dims ,
            theta.k = seq(-5,5,len=15) )
summary(mod1)

#############################################################################
# EXAMPLE 5: Subset of data Example 2
#############################################################################

set.seed(5698)
I <- 11                             # number of items
n <- 3000                           # number of persons
b <- seq(-2,2, len=I)               # item difficulties
theta <- stats::rnorm( n, sd=1.3 )  # person abilities
# define item clusters
itemcluster <- rep(0,I)
itemcluster[ c(3,5)] <- 1
itemcluster[c(2,4,9)] <- 2
# residual correlations
rho <- c( .7 , .5 )
# simulate data
dat <- sirt::sim.rasch.dep( theta , b , itemcluster , rho )
colnames(dat) <- paste("I" , seq(1,ncol(dat)) , sep="")

# select subdataset with only one dependent item cluster
item.sel <- scan( what="character" , nlines=1 )
    I1 I6 I7 I8 I10 I11 I3 I5
dat1 <- dat[,item.sel]

#******************
#*** Model 1a: estimate Copula model in sirt
itemcluster <- rep(0,8)
itemcluster[c(7,8)] <- 1
mod1a <- sirt::rasch.copula2( dat3 , itemcluster=itemcluster )
summary(mod1a)

#******************
#*** Model 1b: estimate Copula model in mirt
library(mirt)
#*** redefine dataset for estimation in mirt
dat2 <- dat1[ , itemcluster == 0 ]
dat2 <- as.data.frame(dat2)
# combine items 3 and 5
dat2$C35 <- dat1[,"I3"] + 2*dat1[,"I5"]
table( dat2$C35 , paste0( dat1[,"I3"],dat1[,"I5"]) )
#* define mirt model
mirtmodel <- mirt::mirt.model("
      F = 1-7
      CONSTRAIN = (1-7,a1)
      " )
#-- Copula function with two dependent items
# define item category function for pseudo-items like C35
P.copula2 <- function(par,Theta, ncat){
     b1 <- par[1]
     b2 <- par[2]
     a1 <- par[3]
     ldelta <- par[4]
     P1 <- stats::plogis( a1*(Theta - b1 ) )
     P2 <- stats::plogis( a1*(Theta - b2 ) )
     Q1 <- 1-P1
     Q2 <- 1-P2
     # define vector-wise minimum function
     minf2 <- function( x1 , x2 ){
         ifelse( x1 < x2 , x1 , x2 )
                                }
     # distribution under independence
     F00 <- Q1*Q2
     F10 <- Q1*Q2 + P1*Q2
     F01 <- Q1*Q2 + Q1*P2
     F11 <- 1+0*Q1
     F_ind <- c(F00,F10,F01,F11)
     # distribution under maximal dependence
     F00 <- minf2(Q1,Q2)
     F10 <- Q2              # = minf2(1,Q2)
     F01 <- Q1              # = minf2(Q1,1)
     F11 <- 1+0*Q1          # = minf2(1,1)
     F_dep <- c(F00,F10,F01,F11)
     # compute mixture distribution
     delta <- stats::plogis(ldelta)
     F_tot <- (1-delta)*F_ind + delta * F_dep
     # recalculate probabilities of mixture distribution
     L1 <- length(Q1)
     v1 <- 1:L1
     F00 <- F_tot[v1]
     F10 <- F_tot[v1+L1]
     F01 <- F_tot[v1+2*L1]
     F11 <- F_tot[v1+3*L1]
     P00 <- F00
     P10 <- F10 - F00
     P01 <- F01 - F00
     P11 <- 1 - F10 - F01 + F00
     prob_tot <- c( P00 , P10 , P01 , P11 )
     return(prob_tot)
        }
# create item
copula2 <- mirt::createItem(name="copula2", par=c(b1 = 0 , b2 = 0.2 , a1=1 , ldelta=0) ,
                est=c(TRUE,TRUE,TRUE,TRUE) , P=P.copula2 ,
                lbound=c(-Inf,-Inf,0,-Inf) , ubound=c(Inf,Inf,Inf,Inf) )
# define item types
itemtype <- c( rep("2PL",6), "copula2" )
customItems <- list("copula2"=copula2)
# parameter table
mod.pars <- mirt::mirt(dat2, 1, itemtype=itemtype,
                customItems=customItems, pars = 'values')
# estimate model
mod1b <- mirt::mirt(dat2, mirtmodel , itemtype=itemtype , customItems=customItems,
                verbose = TRUE , pars=mod.pars ,
                technical=list(customTheta=as.matrix(seq(-4,4,len=21)) ) )
# estimated coefficients
cmod <- sirt::mirt.wrapper.coef(mod)$coef

# compare common item discrimination
round( c("sirt"=mod1a$item$a[1] , "mirt"=cmod$a1[1] ) , 4 )
  ##     sirt   mirt
  ##   1.2845 1.2862
# compare delta parameter
round( c("sirt"=mod1a$item$delta[7] , "mirt"= stats::plogis( cmod$ldelta[7] ) ) , 4 )
  ##     sirt   mirt
  ##   0.6298 0.6297
# compare thresholds a*b
dfr <- cbind( "sirt"=mod1a$item$thresh ,
               "mirt"= c(- cmod$d[-7],cmod$b1[7]*cmod$a1[1]  , cmod$b2[7]*cmod$a1[1]))
round(dfr,4)
  ##           sirt    mirt
  ##   [1,] -1.9236 -1.9231
  ##   [2,] -0.0565 -0.0562
  ##   [3,]  0.3993  0.3996
  ##   [4,]  0.8058  0.8061
  ##   [5,]  1.5293  1.5295
  ##   [6,]  1.9569  1.9572
  ##   [7,] -1.1414 -1.1404
  ##   [8,] -0.4005 -0.3996
}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IRT copula models}
\keyword{Local dependence}
\keyword{summary}
\keyword{logLik}
\keyword{anova}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
