%% File Name: xxirt.Rd
%% File Version: 0.540

\name{xxirt}
\alias{xxirt}
\alias{summary.xxirt}
\alias{print.xxirt}
\alias{logLik.xxirt}
\alias{anova.xxirt}
\alias{coef.xxirt}
\alias{vcov.xxirt}
\alias{confint.xxirt}
\alias{IRT.se.xxirt}
\alias{IRT.expectedCounts.xxirt}
\alias{IRT.irfprob.xxirt}
\alias{IRT.likelihood.xxirt}
\alias{IRT.posterior.xxirt}
\alias{IRT.modelfit.xxirt}
\alias{summary.IRT.modelfit.xxirt}
\alias{IRT.factor.scores.xxirt}
\alias{xxirt_hessian}
\alias{xxirt_sandwich_pml}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
User Defined Item Response Model
}

\description{
Estimates a user defined item response model. Both, item response functions
and latent trait distributions can be specified by the user (see Details).
By default, the EM algorithm is used for estimation. The number of maximum
EM iterations can be defined with the argument \code{maxit}. The \code{xxirt}
function also allows Newton-Raphson optimization by specifying values of maximum
number of iterations in \code{maxit_nr} larger than zero. Typically, a small initial
number of EM iterations should be chosen to obtain reasonable starting values.
}

\usage{
xxirt(dat, Theta=NULL, itemtype=NULL, customItems=NULL, partable=NULL,
       customTheta=NULL, group=NULL, weights=NULL, globconv=1e-06, conv=1e-04,
       maxit=1000, mstep_iter=4, mstep_reltol=1e-06, maxit_nr=0, optimizer_nr="nlminb",
       estimator="ML", control_nr=list(trace=1), h=1E-4, use_grad=TRUE, verbose=TRUE,
       penalty_fun_item=NULL, np_fun_item=NULL, penalty_fun_theta=NULL,
       np_fun_theta=NULL, pml_args=NULL, verbose_index=NULL,
       cv_kfold=0, cv_maxit=10)

\method{summary}{xxirt}(object, digits=3, file=NULL, ...)

\method{print}{xxirt}(x, ...)

\method{anova}{xxirt}(object,...)

\method{coef}{xxirt}(object,...)

\method{logLik}{xxirt}(object,...)

\method{vcov}{xxirt}(object,...)

\method{confint}{xxirt}(object, parm, level=.95, ... )

\method{IRT.expectedCounts}{xxirt}(object,...)

\method{IRT.factor.scores}{xxirt}(object, type="EAP", ...)

\method{IRT.irfprob}{xxirt}(object,...)

\method{IRT.likelihood}{xxirt}(object,...)

\method{IRT.posterior}{xxirt}(object,...)

\method{IRT.modelfit}{xxirt}(object,...)

\method{summary}{IRT.modelfit.xxirt}(object,...)

\method{IRT.se}{xxirt}(object,...)

# computes Hessian matrix
xxirt_hessian(object, h=1e-4, use_shortcut=TRUE)

#- sandwich estimate for pairwise maximum likelihood estimation
xxirt_sandwich_pml(object, h=1e-4)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame with item responses
}
  \item{Theta}{
Matrix with \eqn{\bold{\theta}} grid vector of latent trait
}
  \item{itemtype}{
Vector of item types
}
  \item{customItems}{
List containing types of item response functions created by
\code{\link{xxirt_createDiscItem}}.
}
 \item{partable}{
Item parameter table which is initially created by
\code{\link{xxirt_createParTable}} and which can be modified by
\code{\link{xxirt_modifyParTable}}.
 }
  \item{customTheta}{
User defined \eqn{\bold{\theta}} distribution created by
\code{\link{xxirt_createThetaDistribution}}.
}
  \item{group}{
Optional vector of group indicators
}
  \item{weights}{
Optional vector of person weights
}
  \item{globconv}{
Convergence criterion for relative change in deviance
}
  \item{conv}{
Convergence criterion for absolute change in parameters
}
  \item{maxit}{
Maximum number of iterations in the EM algorithm
}
  \item{mstep_iter}{
Maximum number of iterations in M-step
}
  \item{mstep_reltol}{
Convergence criterion in M-step
}
\item{maxit_nr}{Number of Newton-Raphson iterations after EM algorithm}
\item{optimizer_nr}{Type of optimizer for Newton-Raphson optimization.
Alternatives are \code{"optim"} or \code{"nlminb"} or other options of
\code{\link{sirt_optimizer}}.}
\item{estimator}{Marginal maximum likelihood (\code{"ML"}) or pairwise
maximum likelihood (\code{"PML"})
}
\item{control_nr}{Argument \code{control} for optimizer.}
\item{h}{Numerical differentiation parameter}
\item{use_grad}{Logical indicating whether the gradient should be supplied
to \code{\link[stats:optim]{stats::optim}}}
  \item{verbose}{
Logical indicating whether iteration progress should be displayed
}
\item{penalty_fun_item}{Optional penalty function used in regularized
estimation. Used as a function of \code{x} (vector of item parameters)}
\item{np_fun_item}{Function that counts the number of item parameters in regularized
estimation. Used as a function of \code{x} (vector of item parameters)}
\item{penalty_fun_theta}{Optional penalty function used in regularized
estimation. Used as a function of \code{x} (vector of theta parameters)}
\item{np_fun_theta}{Function that counts the number of theta parameters in regularized
estimation. Used as a function of \code{x} (vector of theta parameters)}
\item{object}{Object of class \code{xxirt}}
\item{digits}{Number of digits to be rounded}
\item{file}{Optional file name to which \code{summary} output is written}
\item{parm}{Optional vector of parameters}
\item{level}{Confidence level}
\item{pml_args}{Weight matrices for \code{estimator="PML"}}
\item{verbose_index}{Logical indicating whether item index should be
printed in estimation output}
\item{cv_kfold}{Number of k folds in cross validation. The default is 0 (no
cross-validation)}
\item{cv_maxit}{Maximum number of iterations for each cross-validation sample}
\item{x}{Object of class \code{xxirt}}
\item{type}{Type of person parameter estimate. Currently, only
\code{EAP} is implemented.}
\item{use_shortcut}{Logical indicating whether a shortcut in the computation
should be utilized}
\item{\dots}{
Further arguments to be passed
}
}

\details{
Item response functions can be specified as functions of unknown parameters
\eqn{\bold{\delta}_i} such that
\eqn{P(X_{i}=x | \bold{\theta})=f_i( x | \bold{\theta} ; \bold{\delta}_i  )}
The item response model is estimated under the assumption of
local stochastic independence of items. Equality constraints of
item parameters \eqn{\bold{\delta}_i} among items are allowed.

The probability distribution \eqn{P(\bold{\theta})} are specified as functions
of an unknown parameter vector \eqn{\bold{\gamma}}.

A penalty function for item parameters can be specified in
\code{penalty_fun_item}. The penalty function should be differentiable and
a non-differentiable function (e.g., the absolute value function) should
be approximated by a differentiable function.
}

\value{
List with following entries

\item{partable}{Item parameter table}
\item{par_items}{Vector with estimated item parameters}
\item{par_items_summary}{Data frame with item parameters}
\item{par_items_bounds}{Data frame with summary on bounds of estimated
 item parameters}
\item{par_Theta}{Vector with estimated parameters of theta distribution}
\item{Theta}{Matrix with \eqn{\bold{\theta}} grid}
\item{probs_items}{Item response functions}
\item{probs_Theta}{Theta distribution}
\item{deviance}{Deviance}
\item{loglik}{Log likelihood value}
\item{ic}{Information criteria}
\item{item_list}{List with item functions}
\item{customItems}{Used customized item response functions}
\item{customTheta}{Used customized theta distribution}
\item{cv_loglike}{Cross-validated log-likelihood value (if \code{cv_kfold>0})}
\item{p.xi.aj}{Individual likelihood}
\item{p.aj.xi}{Individual posterior}
\item{ll_case}{Case-wise log-likelihood values}
\item{n.ik}{Array of expected counts}
\item{EAP}{EAP person parameter estimates}
\item{dat}{Used dataset with item responses}
\item{dat_resp}{Dataset with response indicators}
\item{weights}{Vector of person weights}
\item{G}{Number of groups}
\item{group}{Integer vector of group indicators}
\item{group_orig}{Vector of original group_identifiers}
\item{ncat}{Number of categories per item}
\item{converged}{Logical whether model has converged}
\item{iter}{Number of iterations needed}
}

%\references{
%% ~put references to the literature/web site here ~
%}



%\author{
%Alexander Robitzsch
%}



%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See the \code{\link[mirt:createItem]{mirt::createItem}} and
\code{\link[mirt:mirt]{mirt::mirt}} functions in the \pkg{mirt}
package for similar functionality.
}

\examples{
\dontrun{
#############################################################################
## EXAMPLE 1: Unidimensional item response functions
#############################################################################

data(data.read)
dat <- data.read

#------ Definition of item response functions

#*** IRF 2PL
P_2PL <- function( par, Theta, ncat){
    a <- par[1]
    b <- par[2]
    TP <- nrow(Theta)
    P <- matrix( NA, nrow=TP, ncol=ncat)
    P[,1] <- 1
    for (cc in 2:ncat){
        P[,cc] <- exp( (cc-1) * a * Theta[,1] - b )
    }
    P <- P / rowSums(P)
    return(P)
}

#*** IRF 1PL
P_1PL <- function( par, Theta, ncat){
    b <- par[1]
    TP <- nrow(Theta)
    P <- matrix( NA, nrow=TP, ncol=ncat)
    P[,1] <- 1
    for (cc in 2:ncat){
        P[,cc] <- exp( (cc-1) * Theta[,1] - b )
    }
    P <- P / rowSums(P)
    return(P)
}

#** created item classes of 1PL and 2PL models
par <- c( "a"=1, "b"=0 )
# define some slightly informative prior of 2PL
item_2PL <- sirt::xxirt_createDiscItem( name="2PL", par=par, est=c(TRUE,TRUE),
               P=P_2PL, prior=c(a="dlnorm"), prior_par1=c( a=0 ),
               prior_par2=c(a=5) )
item_1PL <- sirt::xxirt_createDiscItem( name="1PL", par=par[2], est=c(TRUE),
               P=P_1PL )
customItems <- list( item_1PL,  item_2PL )

#---- definition theta distribution

#** theta grid
Theta <- matrix( seq(-6,6,length=21), ncol=1 )

#** theta distribution
P_Theta1 <- function( par, Theta, G){
    mu <- par[1]
    sigma <- max( par[2], .01 )
    TP <- nrow(Theta)
    pi_Theta <- matrix( 0, nrow=TP, ncol=G)
    pi1 <- dnorm( Theta[,1], mean=mu, sd=sigma )
    pi1 <- pi1 / sum(pi1)
    pi_Theta[,1] <- pi1
    return(pi_Theta)
}
#** create distribution class
par_Theta <- c( "mu"=0, "sigma"=1 )
customTheta  <- sirt::xxirt_createThetaDistribution( par=par_Theta, est=c(FALSE,TRUE),
                       P=P_Theta1 )

#****************************************************************************
#******* Model 1: Rasch model

#-- create parameter table
itemtype <- rep( "1PL", 12 )
partable <- sirt::xxirt_createParTable( dat, itemtype=itemtype,
                        customItems=customItems )

# estimate model
mod1 <- sirt::xxirt( dat=dat, Theta=Theta, partable=partable,
                   customItems=customItems, customTheta=customTheta)
summary(mod1)

# estimate Rasch model by providing starting values
partable1 <- sirt::xxirt_modifyParTable( partable, parname="b",
                   value=- stats::qlogis( colMeans(dat) ) )
# estimate model again
mod1b <- sirt::xxirt( dat=dat, Theta=Theta, partable=partable1,
                   customItems=customItems, customTheta=customTheta )
summary(mod1b)

# extract coefficients, covariance matrix and standard errors
coef(mod1b)
vcov(mod1b)
IRT.se(mod1b)

#** start with EM and finalize with Newton-Raphson algorithm
mod1c <- sirt::xxirt( dat=dat, Theta=Theta, partable=partable,
                   customItems=customItems, customTheta=customTheta,
                   maxit=20, maxit_nr=300)
summary(mod1c)

#****************************************************************************
#******* Model 2: 2PL Model with three groups of item discriminations

#-- create parameter table
itemtype <- rep( "2PL", 12 )
partable <- sirt::xxirt_createParTable( dat, itemtype=itemtype, customItems=customItems)
# modify parameter table: set constraints for item groups A, B and C
partable1 <- sirt::xxirt_modifyParTable(partable, item=paste0("A",1:4),
                         parname="a", parindex=111)
partable1 <- sirt::xxirt_modifyParTable(partable1, item=paste0("B",1:4),
                         parname="a", parindex=112)
partable1 <- sirt::xxirt_modifyParTable(partable1, item=paste0("C",1:4),
                         parname="a", parindex=113)
# delete prior distributions
partable1 <- sirt::xxirt_modifyParTable(partable1, parname="a", prior=NA)

#-- fix sigma to 1
customTheta1 <- customTheta
customTheta1$est <- c("mu"=FALSE,"sigma"=FALSE )

# estimate model
mod2 <- sirt::xxirt( dat=dat, Theta=Theta, partable=partable1,
                  customItems=customItems, customTheta=customTheta1 )
summary(mod2)

#****************************************************************************
#******* Model 3: Cloglog link function

#*** IRF cloglog
P_1N <- function( par, Theta, ncat){
    b <- par
    TP <- nrow(Theta)
    P <- matrix( NA, nrow=TP, ncol=ncat)
    P[,2] <- 1 - exp( - exp( Theta - b ) )
    P[,1] <- 1 - P[,2]
    return(P)
}
par <- c("b"=0)
item_1N <- sirt::xxirt_createDiscItem( name="1N", par=par, est=c(TRUE),
                    P=P_1N )
customItems <- list( item_1N )
itemtype <- rep( "1N", I )
partable <- sirt::xxirt_createParTable( dat[,items], itemtype=itemtype,
                      customItems=customItems )
partable <- sirt::xxirt_modifyParTable( partable=partable, parname="b",
                 value=- stats::qnorm( colMeans(dat[,items] )) )

#*** estimate model
mod3 <- sirt::xxirt( dat=dat, Theta=Theta, partable=partable, customItems=customItems,
                customTheta=customTheta )
summary(mod3)
IRT.compareModels(mod1,mod3)

#****************************************************************************
#******* Model 4: Latent class model

K <- 3 # number of classes
Theta <- diag(K)

#*** Theta distribution
P_Theta1 <- function( par, Theta, G  ){
    logitprobs <- par[1:(K-1)]
    l1 <- exp( c( logitprobs, 0 ) )
    probs <- matrix( l1/sum(l1), ncol=1)
    return(probs)
}

par_Theta <- stats::qlogis( rep( 1/K, K-1 ) )
names(par_Theta) <- paste0("pi",1:(K-1) )
customTheta  <- sirt::xxirt_createThetaDistribution( par=par_Theta,
                     est=rep(TRUE,K-1), P=P_Theta1)

#*** IRF latent class
P_lc <- function( par, Theta, ncat){
    b <- par
    TP <- nrow(Theta)
    P <- matrix( NA, nrow=TP, ncol=ncat)
    P[,1] <- 1
    for (cc in 2:ncat){
        P[,cc] <- exp( Theta \%*\% b )
    }
    P <- P / rowSums(P)
    return(P)
}
par <- seq( -1.5, 1.5, length=K )
names(par) <- paste0("b",1:K)
item_lc <- sirt::xxirt_createDiscItem( name="LC", par=par,
                 est=rep(TRUE,K), P=P_lc )
customItems <- list( item_lc )

# create parameter table
itemtype <- rep( "LC", 12 )
partable <- sirt::xxirt_createParTable( dat, itemtype=itemtype, customItems=customItems)
partable

#*** estimate model
mod4 <- sirt::xxirt( dat=dat, Theta=Theta, partable=partable, customItems=customItems,
                customTheta=customTheta)
summary(mod4)
# class probabilities
mod4$probs_Theta
# item response functions
imod4 <- IRT.irfprob( mod5 )
round( imod4[,2,], 3 )

#****************************************************************************
#******* Model 5: Ordered latent class model

K <- 3 # number of classes
Theta <- diag(K)
Theta <- apply( Theta, 1, cumsum )

#*** Theta distribution
P_Theta1 <- function( par, Theta, G  ){
    logitprobs <- par[1:(K-1)]
    l1 <- exp( c( logitprobs, 0 ) )
    probs <- matrix( l1/sum(l1), ncol=1)
    return(probs)
}
par_Theta <- stats::qlogis( rep( 1/K, K-1 ) )
names(par_Theta) <- paste0("pi",1:(K-1) )
customTheta  <- sirt::xxirt_createThetaDistribution( par=par_Theta,
                est=rep(TRUE,K-1), P=P_Theta1  )

#*** IRF ordered latent class
P_olc <- function( par, Theta, ncat){
    b <- par
    TP <- nrow(Theta)
    P <- matrix( NA, nrow=TP, ncol=ncat)
    P[,1] <- 1
    for (cc in 2:ncat){
        P[,cc] <- exp( Theta \%*\% b )
    }
    P <- P / rowSums(P)
    return(P)
}

par <- c( -1, rep( .5,, length=K-1 ) )
names(par) <- paste0("b",1:K)
item_olc <- sirt::xxirt_createDiscItem( name="OLC", par=par, est=rep(TRUE,K),
                    P=P_olc, lower=c( -Inf, 0, 0 ) )
customItems <- list( item_olc )
itemtype <- rep( "OLC", 12 )
partable <- sirt::xxirt_createParTable( dat, itemtype=itemtype, customItems=customItems)
partable

#*** estimate model
mod5 <- sirt::xxirt( dat=dat, Theta=Theta, partable=partable, customItems=customItems,
                customTheta=customTheta )
summary(mod5)
# estimated item response functions
imod5 <- IRT.irfprob( mod5 )
round( imod5[,2,], 3 )

#############################################################################
## EXAMPLE 2: Multiple group models with xxirt
#############################################################################

data(data.math)
dat <- data.math$data
items <- grep( "M[A-Z]", colnames(dat), value=TRUE )
I <- length(items)

Theta <- matrix( seq(-8,8,len=31), ncol=1 )

#****************************************************************************
#******* Model 1: Rasch model, single group

#*** Theta distribution
P_Theta1 <- function( par, Theta, G  ){
    mu <- par[1]
    sigma <- max( par[2], .01 )
    p1 <- stats::dnorm( Theta[,1], mean=mu, sd=sigma)
    p1 <- p1 / sum(p1)
    probs <- matrix( p1, ncol=1)
    return(probs)
}

par_Theta <- c(0,1)
names(par_Theta) <- c("mu","sigma")
customTheta  <- sirt::xxirt_createThetaDistribution( par=par_Theta,
                   est=c(FALSE,TRUE), P=P_Theta1  )
customTheta

#*** IRF 1PL logit
P_1PL <- function( par, Theta, ncat){
    b <- par
    TP <- nrow(Theta)
    P <- matrix( NA, nrow=TP, ncol=ncat)
    P[,2] <- plogis( Theta - b )
    P[,1] <- 1 - P[,2]
    return(P)
}
par <- c("b"=0)
item_1PL <- sirt::xxirt_createDiscItem( name="1PL", par=par, est=c(TRUE), P=P_1PL)
customItems <- list( item_1PL )

itemtype <- rep( "1PL", I )
partable <- sirt::xxirt_createParTable( dat[,items], itemtype=itemtype,
                       customItems=customItems )
partable <- sirt::xxirt_modifyParTable( partable=partable, parname="b",
                  value=- stats::qlogis( colMeans(dat[,items] )) )

#*** estimate model
mod1 <- sirt::xxirt( dat=dat[,items], Theta=Theta, partable=partable,
                customItems=customItems, customTheta=customTheta )
summary(mod1)

#****************************************************************************
#******* Model 2: Rasch model, multiple groups

#*** Theta distribution
P_Theta2 <- function( par, Theta, G  ){
    mu1 <- par[1]
    mu2 <- par[2]
    sigma1 <- max( par[3], .01 )
    sigma2 <- max( par[4], .01 )
    TP <- nrow(Theta)
    probs <- matrix( NA, nrow=TP, ncol=G)
    p1 <- stats::dnorm( Theta[,1], mean=mu1, sd=sigma1)
    probs[,1] <- p1 / sum(p1)
    p1 <- stats::dnorm( Theta[,1], mean=mu2, sd=sigma2)
    probs[,2] <- p1 / sum(p1)
    return(probs)
}
par_Theta <- c(0,0,1,1)
names(par_Theta) <- c("mu1","mu2","sigma1","sigma2")
customTheta2  <- sirt::xxirt_createThetaDistribution( par=par_Theta,
                    est=c(FALSE,TRUE,TRUE,TRUE), P=P_Theta2  )
print(customTheta2)

#*** estimate model
mod2 <- sirt::xxirt( dat=dat[,items], group=dat$female, Theta=Theta, partable=partable,
           customItems=customItems, customTheta=customTheta2, maxit=40)
summary(mod2)
IRT.compareModels(mod1, mod2)

#*** compare results with TAM package
library(TAM)
mod2b <- TAM::tam.mml( resp=dat[,items], group=dat$female )
summary(mod2b)
IRT.compareModels(mod1, mod2, mod2b)

#############################################################################
## EXAMPLE 3: Regularized 2PL model
#############################################################################

data(data.read, package="sirt")
dat <- data.read

#------ Definition of item response functions

#*** IRF 2PL
P_2PL <- function( par, Theta, ncat){
    a <- par[1]
    b <- par[2]
    TP <- nrow(Theta)
    P <- matrix( NA, nrow=TP, ncol=ncat)
    P[,1] <- 1
    for (cc in 2:ncat){
        P[,cc] <- exp( (cc-1) * a * Theta[,1] - b )
    }
    P <- P / rowSums(P)
    return(P)
}

#** created item classes of 1PL and 2PL models
par <- c( "a"=1, "b"=0 )
# define some slightly informative prior of 2PL
item_2PL <- sirt::xxirt_createDiscItem( name="2PL", par=par, est=c(TRUE,TRUE),
               P=P_2PL, prior=c(a="dlnorm"), prior_par1=c( a=0 ),
               prior_par2=c(a=5) )
customItems <- list( item_2PL )

#---- definition theta distribution

#** theta grid
Theta <- matrix( seq(-6,6,length=21), ncol=1 )

#** theta distribution
P_Theta1 <- function( par, Theta, G){
    mu <- par[1]
    sigma <- max( par[2], .01 )
    TP <- nrow(Theta)
    pi_Theta <- matrix( 0, nrow=TP, ncol=G)
    pi1 <- dnorm( Theta[,1], mean=mu, sd=sigma )
    pi1 <- pi1 / sum(pi1)
    pi_Theta[,1] <- pi1
    return(pi_Theta)
}
#** create distribution class
par_Theta <- c( "mu"=0, "sigma"=1 )
customTheta  <- sirt::xxirt_createThetaDistribution( par=par_Theta, est=c(FALSE,FALSE),
                       P=P_Theta1 )

#****************************************************************************
#******* Model 1: 2PL model

itemtype <- rep( "2PL", 12 )
partable <- sirt::xxirt_createParTable( dat, itemtype=itemtype,
                        customItems=customItems )

mod1 <- sirt::xxirt( dat=dat, Theta=Theta, partable=partable,
                   customItems=customItems, customTheta=customTheta)
summary(mod1)

#****************************************************************************
#******* Model 2: Regularized 2PL model with regularization on item loadings

# define regularized estimation of item loadings
parindex <- partable[ partable$parname=="a","parindex"]

#** penalty is defined by -N*lambda*sum_i (a_i-1)^2
N <- nrow(dat)
lambda <- .02
penalty_fun_item <- function(x)
{
    val <- N*lambda*sum( ( x[parindex]-1)^2)
    return(val)
}
# estimate standard deviation
customTheta1  <- sirt::xxirt_createThetaDistribution( par=par_Theta, est=c(FALSE,TRUE),
                       P=P_Theta1 )
mod2 <- sirt::xxirt( dat=dat, Theta=Theta, partable=partable,
                   customItems=customItems, customTheta=customTheta1,
                   penalty_fun_item=penalty_fun_item)
summary(mod2)

#############################################################################
## EXAMPLE 4: 2PL mixture model
#############################################################################

#*** simulate data
set.seed(123)
N <- 4000   # number of persons
I <- 15     # number of items
prop <- .25 # mixture proportion for second class

# discriminations and difficulties in first class
a1 <- rep(1,I)
b1 <- seq(-2,2,len=I)
# distribution in second class
mu2 <- 1
sigma2 <- 1.2
# compute parameters with constraint N(0,1) in second class
# a*(sigma*theta+mu-b)=a*sigma*(theta-(b-mu)/sigma)
#=> a2=a*sigma and b2=(b-mu)/sigma
a2 <- a1
a2[c(2,4,6,8)] <- 0.2  # some items with different discriminations
a2 <- a2*sigma2
b2 <- b1
b2[1:5] <- 1   # first 5 item with different difficulties
b2 <- (b2-mu2)/sigma2
dat1 <- sirt::sim.raschtype(theta=stats::rnorm(N*(1-prop)), b=b1, fixed.a=a1)
dat2 <- sirt::sim.raschtype(theta=stats::rnorm(N*prop), b=b2, fixed.a=a2)
dat <- rbind(dat1, dat2)

#**** model specification

#*** define theta distribution
TP <- 21
theta <- seq(-6,6,length=TP)
# stack theta vectors below each others=> 2 latent classes
Theta <- matrix( c(theta, theta ), ncol=1 )
# distribution of theta (i.e., N(0,1))
w_theta <- dnorm(theta)
w_theta <- w_theta / sum(w_theta)

P_Theta1 <- function( par, Theta, G){
    p2_logis <- par[1]
    p2 <- stats::plogis( p2_logis )
    p1 <- 1-p2
    pi_Theta <- c( p1*w_theta, p2*w_theta)
    pi_Theta <- matrix(pi_Theta, ncol=1)
    return(pi_Theta)
}

par_Theta <- c( p2_logis=qlogis(.25))
customTheta  <- sirt::xxirt_createThetaDistribution( par=par_Theta, est=c(TRUE),
                       P=P_Theta1)

# IRF for 2-class mixture 2PL model
par <- c(a1=1, a2=1, b1=0, b2=.5)

P_2PLmix <- function( par, Theta, ncat)
{
    a1 <- par[1]
    a2 <- par[2]
    b1 <- par[3]
    b2 <- par[4]
    P <- matrix( NA, nrow=2*TP, ncol=ncat)
    TP <- nrow(Theta)/2
    P1 <- stats::plogis( a1*(Theta[1:TP,1]-b1) )
    P2 <- stats::plogis( a2*(Theta[TP+1:(2*TP),1]-b2) )
    P[,2] <- c(P1, P2)
    P[,1] <- 1-P[,2]
    return(P)
}

# define some slightly informative prior of 2PL
item_2PLmix <- sirt::xxirt_createDiscItem( name="2PLmix", par=par,
               est=c(TRUE,TRUE,TRUE,TRUE), P=P_2PLmix )
customItems <- list( item_2PLmix )

#****************************************************************************
#******* Model 1: 2PL mixture model

itemtype <- rep( "2PLmix", I )
partable <- sirt::xxirt_createParTable( dat, itemtype=itemtype,
                        customItems=customItems )
mod1 <- sirt::xxirt( dat=dat, Theta=Theta, partable=partable,
                   customItems=customItems, customTheta=customTheta)
summary(mod1)

#############################################################################
## EXAMPLE 5: Partial credit model with MML and PML
#############################################################################

data(data.gpcm, package="TAM")
dat <- data.gpcm

# recode data and include some missings
dat[ dat[,1]==3, 1] <- 2
dat[ 1, c(1,2) ] <- NA
dat[2,3] <- NA

#------ Definition of item response functions
#*** IRF 2PL
P_2PL <- function( par, Theta, ncat)
{
    b <- par
    TP <- nrow(Theta)
    P <- matrix( NA, nrow=TP, ncol=ncat)
    P[,1] <- 1
    for (cc in 2:ncat){
        P[,cc] <- exp( Theta[,1] - b[cc-1] )
    }
    P <- P / rowSums(P)
    return(P)
}

P_PCM2 <- function( par, Theta, ncat=3)
{
    P <- P_2PL(par=par, Theta=Theta, ncat=ncat)
    return(P)
}

P_PCM3 <- function( par, Theta, ncat=4)
{
    P <- P_2PL(par=par, Theta=Theta, ncat=ncat)
    return(P)
}

# define some slightly informative prior of 2PL
par <- c( b1=-1, b2=1 )
NP <- length(par)
item_PCM2 <- sirt::xxirt_createDiscItem( name="PCM2", par=par, est=rep(TRUE,NP),
               P=P_PCM2 )
par <- c( b1=-1, b2=0, b3=1 ) ; NP <- length(par)
item_PCM3 <- sirt::xxirt_createDiscItem( name="PCM3", par=par, est=rep(TRUE,NP),
               P=P_PCM3 )

customItems <- list( item_PCM2,  item_PCM3 )

#---- definition theta distribution
#** theta grid
Theta <- matrix( seq(-6,6,length=21), ncol=1 )

#** theta distribution
P_Theta1 <- function( par, Theta, G){
    mu <- par[1]
    sigma <- max( par[2], .01 )
    TP <- nrow(Theta)
    pi_Theta <- matrix( 0, nrow=TP, ncol=G)
    pi1 <- dnorm( Theta[,1], mean=mu, sd=sigma )
    pi1 <- pi1 / sum(pi1)
    pi_Theta[,1] <- pi1
    return(pi_Theta)
}
#** create distribution class
par_Theta <- c( "mu"=0, "sigma"=1 )
customTheta  <- sirt::xxirt_createThetaDistribution( par=par_Theta, est=c(FALSE,TRUE),
                       P=P_Theta1 )

#-- create parameter table
itemtype <- c("PCM2", "PCM3", "PCM3")
partable <- sirt::xxirt_createParTable( dat, itemtype=itemtype,
                        customItems=customItems )

#***** Model 1: MML
mod1 <- sirt::xxirt( dat=dat, Theta=Theta, partable=partable,
                   customItems=customItems, customTheta=customTheta)
summary(mod1)

#***** Model 2: PML
I <- ncol(dat)
W1 <- rep(1,I)
W2 <- 1-diag(I)
W2[3,2] <- 0
W2[ upper.tri(W2)] <- 0
pml_args <- list(W1=W1/sum(W1), W2=W2/sum(W2) )

mod2 <- sirt::xxirt( dat=dat, Theta=Theta, partable=partable,
                   customItems=customItems, customTheta=mod1$customTheta,
                   estimator="PML", pml_args=pml_args)

# variance matrix
smod2 <- sirt::xxirt_sandwich_pml(object=mod2)
smod2


#############################################################################
## EXAMPLE 6: Person covariates for item parameters and trait distribution
#############################################################################

#-- simulate data
set.seed(2025)
N <- 2000
I <- 10
b <- c(0, seq(-2,2,len=I-1))
b1 <- b
b1[1] <- 1
dat1 <- sirt::sim.raschtype( stats::rnorm(N, mean=0, sd=1), b=b )
dat2 <- sirt::sim.raschtype( stats::rnorm(N, mean=1, sd=1.5), b=b1 )
dat <- data.frame( group=rep( c(1,2), each=N), rbind(dat1, dat2))

#***************************************************************************
#***** Model 1: two-group Rasch model

#*** IRF 1PL
P_1PL <- function( par, Theta, ncat){
    b <- par[1]
    TP <- nrow(Theta)
    P <- matrix( NA, nrow=TP, ncol=ncat)
    P[,1] <- 1
    for (cc in 2:ncat){
        P[,cc] <- exp( (cc-1) * Theta[,1] - b )
    }
    P <- P / rowSums(P)
    return(P)
}

#** create item classes of 1PL and 2PL models
par <- c( b=0 )
item_1PL <- sirt::xxirt_createDiscItem( name="1PL", par=par, est=c(TRUE), P=P_1PL )
customItems <- list( item_1PL )

#---- definition theta distribution

#** theta grid
Theta <- matrix( seq(-6,6,length=21), ncol=1 )

#** theta distribution
P_Theta1 <- function( par, Theta, G)
{
    mu1 <- par[1]
    sigma1 <- max( par[2], .01 )
    mu2 <- par[3]
    sigma2 <- max( par[4], .01 )
    TP <- nrow(Theta)
    pi_Theta <- matrix( 0, nrow=TP, ncol=G)
    pi1 <- dnorm( Theta[,1], mean=mu1, sd=sigma1 )
    pi1 <- pi1 / sum(pi1)
    pi_Theta[,1] <- pi1
    pi1 <- dnorm( Theta[,1], mean=mu2, sd=sigma2 )
    pi1 <- pi1 / sum(pi1)
    pi_Theta[,2] <- pi1
    return(pi_Theta)
}
#** create distribution class
par_Theta <- c( "mu1"=0, "sigma1"=1, "mu2"=0, "sigma2"=1 )
customTheta  <- sirt::xxirt_createThetaDistribution( par=par_Theta,
                       est=c(FALSE,TRUE,TRUE,TRUE), P=P_Theta1 )

#-- create parameter table
itemtype <- rep( "1PL", I )
partable <- sirt::xxirt_createParTable( dat[,-1], itemtype=itemtype,
                        customItems=customItems )

# estimate model
mod1 <- sirt::xxirt( dat=dat[,-1], Theta=Theta, partable=partable, group=dat$group,
                   customItems=customItems, customTheta=customTheta)
summary(mod1)

#***************************************************************************
#***** Model 2: Rasch model with person covariates in trait distribution

#- matrix with person covariates
X <- dat[,"group", drop=FALSE]

#** theta distribution
P_Theta2 <- function( par, Theta, G, X)
{
    mu1 <- par[1]
    sigma1 <- max( par[2], .01 )
    mu2 <- par[3]
    sigma2 <- max( par[4], .01 )
    TP <- nrow(Theta)
    N <- nrow(X)
    pi_Theta <- matrix( 0, nrow=TP, ncol=N)
    mu <- par[ 2*(X[,1]-1)+1 ]
    sigma <- par[ 2*(X[,1]-1)+2 ]
    for (tt in 1:TP){
        pi_Theta[tt,] <- stats::dnorm(Theta[tt,1], mean=mu, sd=sigma)
    }
    cp <- sirt::sirt_matrix2( colSums(pi_Theta), nrow=TP)
    pi_Theta <- pi_Theta / ( cp + 1e-80 )
    return(pi_Theta)
}

par_Theta2 <- mod1$customTheta$par

customTheta2  <- sirt::xxirt_createThetaDistribution( par=par_Theta2,
                       est=c(FALSE,TRUE,TRUE,TRUE),
                       P=P_Theta2, X=X, lower=c(-10,0.01, -10, 0.01) )

partable2 <- mod1$partable

mod2 <- sirt::xxirt( dat=dat[,-1], Theta=Theta, partable=partable2,
                   customItems=customItems, customTheta=customTheta2)
summary(mod2)

#***************************************************************************
#***** Model 3: Rasch model with person covariates in trait distribution and
#*****          item parameters

#*** IRF 1PL
P_1PL <- function( par, Theta, ncat, X)
{
    TP <- nrow(Theta)
    N <- nrow(X)
    P <- array( 0, dim=c(TP,ncat,N) )
    efficient <- TRUE
    if (! efficient ){
        bM <- sirt::sirt_matrix2( par[ X[,1] ], nrow=TP)
        P1 <- stats::plogis(Theta[,1] - bM )
    } else {
        h1 <- exp(-Theta[,1])
        h2 <- sirt::sirt_matrix2( exp(par[ X[,1] ]), nrow=TP)
        P1 <- 1/ ( 1+h1*h2 )
    }
    P[,2,] <- P1
    P[,1,] <- 1 - P[,2,]
    return(P)
}

#** create item classes of 1PL and 2PL models
par <- c( b1=0, b2=-1 )
item_1PL <- sirt::xxirt_createDiscItem( name="1PL", par=par, est=c(TRUE,TRUE),
                            P=P_1PL, X=X )

customItems <- list( item_1PL )

#-- create parameter table
itemtype <- rep( "1PL", I )
partable2 <- sirt::xxirt_createParTable( dat[,-1], itemtype=itemtype,
                        customItems=customItems )

#-- distribution parameters
customTheta2b <- mod2$customTheta
customTheta2b$par[3] <- 0      # fix mu2 to 0
customTheta2b$est <- c(FALSE,TRUE,FALSE,TRUE)

# estimate model
mod3 <- sirt::xxirt( dat=resp, Theta=Theta, partable=partable2,
                   customItems=customItems, customTheta=customTheta2b, mstep_iter=2)
summary(mod3)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
