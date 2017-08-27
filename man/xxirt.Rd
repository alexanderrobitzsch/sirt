%% File Name: xxirt.Rd
%% File Version: 0.41
%% File Last Change: 2017-03-01 19:30:55

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

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
User Defined Item Response Model
}

\description{
Estimates a user defined item response model. Both, item response functions
and latent trait distributions can be specified by the user (see Details).
}

\usage{
xxirt(dat, Theta = NULL, itemtype = NULL, customItems = NULL, partable = NULL, 
       customTheta = NULL, group = NULL, weights = NULL, globconv = 1e-06, conv = 1e-04,
       maxit = 200, mstep_iter = 4, mstep_reltol = 1e-06, h = 1E-4 , use_grad = TRUE ,
       verbose = TRUE)
    
\method{summary}{xxirt}(object, digits = 3, file = NULL, ...)

\method{print}{xxirt}(x, ...)

\method{anova}{xxirt}(object,...)

\method{coef}{xxirt}(object,...)

\method{logLik}{xxirt}(object,...)

\method{vcov}{xxirt}(object,...)

\method{confint}{xxirt}(object, parm, level = .95, ... )

\method{IRT.expectedCounts}{xxirt}(object,...)

\method{IRT.factor.scores}{xxirt}(object, type = "EAP", ...)

\method{IRT.irfprob}{xxirt}(object,...)

\method{IRT.likelihood}{xxirt}(object,...)

\method{IRT.posterior}{xxirt}(object,...)

\method{IRT.modelfit}{xxirt}(object,...)

\method{summary}{IRT.modelfit.xxirt}(object,...)

\method{IRT.se}{xxirt}(object,...)

# computes Hessian matrix
xxirt_hessian(object)
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
Maximum number of iterations
}
  \item{mstep_iter}{
Maximum number of iterations in M-step
}
  \item{mstep_reltol}{
Convergence criterion in M-step
}
\item{h}{Numerical differentiation parameter}
\item{use_grad}{Logical indicating whether the gradient should be supplied
to \code{\link[stats:optim]{stats::optim}}}
  \item{verbose}{
Logical indicating whether iteration progress should be displayed
}
\item{object}{Object of class \code{xxirt}}
\item{digits}{Number of digits to be rounded}
\item{file}{Optional file name to which \code{summary} output is written}
\item{parm}{Optional vector of parameters}
\item{level}{Confidence level}
\item{x}{Object of class \code{xxirt}}
\item{type}{Type of person parameter estimate. Currently, only
\code{EAP} is implemented.}
\item{\dots}{
Further arguments to be passed
}
}

\details{
Item response functions can be specified as functions of unknown parameters 
\eqn{\bold{\delta}_i} such that 
\eqn{P(X_{i}=x | \bold{\theta}) = f_i( x | \bold{\theta} ; \bold{\delta}_i  )}
The item response model is estimated under the assumption of
local stochastic independence of items. Equality constraints of
item parameters \eqn{\bold{\delta}_i} among items are allowed. 

Probability distribution \eqn{P(\bold{\theta})} are specified as functions
of an unknown parameter vector \eqn{\bold{\gamma}}.
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
\item{p.xi.aj}{Individual likelihood}
\item{p.aj.xi}{Individual posterior}
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


\author{
Alexander Robitzsch
}

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
P_2PL <- function( par, Theta , ncat){
    a <- par[1]
    b <- par[2]    
    TP <- nrow(Theta)
    P <- matrix( NA , nrow=TP , ncol=ncat)
    P[,1] <- 1
    for (cc in 2:ncat){
        P[,cc] <- exp( (cc-1) * a * Theta[,1] - b )
    }
    P <- P / rowSums(P)    
    return(P)
}

#*** IRF 1PL
P_1PL <- function( par, Theta , ncat){
    b <- par[1]    
    TP <- nrow(Theta)
    P <- matrix( NA , nrow=TP , ncol=ncat)
    P[,1] <- 1
    for (cc in 2:ncat){
        P[,cc] <- exp( (cc-1) * Theta[,1] - b )
    }
    P <- P / rowSums(P)    
    return(P)
}

#** created item classes of 1PL and 2PL models
par <- c( "a"= 1 , "b" = 0 )
# define some slightly informative prior of 2PL
item_2PL <- sirt::xxirt_createDiscItem( name = "2PL" , par = par , est = c(TRUE,TRUE), 
               P = P_2PL, prior = c(a = "dlnorm"), prior_par1 = c( a = 0 ), 
               prior_par2 = c(a=5) )
item_1PL <- sirt::xxirt_createDiscItem( name = "1PL" , par = par[2] , est = c(TRUE), 
               P = P_1PL )
customItems <- list( item_1PL ,  item_2PL )

#---- definition theta distribution

#** theta grid
Theta <- matrix( seq(-6,6,length=21) , ncol=1 )

#** theta distribution
P_Theta1 <- function( par , Theta , G){
    mu <- par[1]
    sigma <- max( par[2] , .01 )
    TP <- nrow(Theta)
    pi_Theta <- matrix( 0 , nrow=TP , ncol=G)
    pi1 <- dnorm( Theta[,1] , mean = mu , sd = sigma )
    pi1 <- pi1 / sum(pi1)
    pi_Theta[,1] <- pi1
    return(pi_Theta)
}
#** create distribution class
par_Theta <- c( "mu"=0, "sigma" = 1 )  
customTheta  <- sirt::xxirt_createThetaDistribution( par=par_Theta , est=c(FALSE,TRUE), 
                       P=P_Theta1 )

#****************************************************************************
#******* Model 1: Rasch model

#-- create parameter table
itemtype <- rep( "1PL" , 12 )
partable <- sirt::xxirt_createParTable( dat , itemtype = itemtype, 
                        customItems = customItems )

# estimate model
mod1 <- sirt::xxirt( dat = dat , Theta=Theta , partable = partable, 
                   customItems = customItems, customTheta = customTheta)
summary(mod1)

# estimate Rasch model by providing starting values
partable1 <- sirt::xxirt_modifyParTable( partable , parname = "b", 
                   value = - stats::qlogis( colMeans(dat) ) )
# estimate model again
mod1b <- sirt::xxirt( dat = dat , Theta=Theta , partable = partable1, 
                   customItems = customItems, customTheta = customTheta )
summary(mod1b)

# extract coefficients, covariance matrix and standard errors
coef(mod1b)
vcov(mod1b)
IRT.se(mod1b)

#****************************************************************************
#******* Model 2: 2PL Model with three groups of item discriminations

#-- create parameter table
itemtype <- rep( "2PL" , 12 )
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
mod2 <- sirt::xxirt( dat = dat, Theta=Theta, partable = partable1, 
                  customItems=customItems, customTheta = customTheta1 )
summary(mod2)

#****************************************************************************
#******* Model 3: Cloglog link function

#*** IRF cloglog
P_1N <- function( par, Theta , ncat){
    b <- par
    TP <- nrow(Theta)
    P <- matrix( NA , nrow=TP , ncol=ncat)
    P[,2] <- 1 - exp( - exp( Theta - b ) )
    P[,1] <- 1 - P[,2]
    return(P)
}
par <- c("b"=0)
item_1N <- sirt::xxirt_createDiscItem( name = "1N" , par = par , est = c(TRUE), 
                    P = P_1N )
customItems <- list( item_1N )
itemtype <- rep( "1N" , I )
partable <- sirt::xxirt_createParTable( dat[,items] , itemtype = itemtype, 
                      customItems = customItems )
partable <- sirt::xxirt_modifyParTable( partable=partable , parname = "b" , 
                 value = - stats::qnorm( colMeans(dat[,items] )) )

#*** estimate model
mod3 <- sirt::xxirt( dat=dat, Theta=Theta, partable=partable, customItems=customItems,
                customTheta= customTheta )                               
summary(mod3)    
IRT.compareModels(mod1,mod3)

#****************************************************************************
#******* Model 4: Latent class model

K <- 3 # number of classes
Theta <- diag(K)

#*** Theta distribution
P_Theta1 <- function( par , Theta , G  ){
    logitprobs <- par[1:(K-1)]
    l1 <- exp( c( logitprobs , 0 ) )    
    probs <- matrix( l1/sum(l1) , ncol=1)    
    return(probs)
}
               
par_Theta <- stats::qlogis( rep( 1/K , K-1 ) )
names(par_Theta) <- paste0("pi",1:(K-1) )
customTheta  <- sirt::xxirt_createThetaDistribution( par=par_Theta, 
                     est= rep(TRUE,K-1) , P=P_Theta1)

#*** IRF latent class
P_lc <- function( par, Theta , ncat){
    b <- par
    TP <- nrow(Theta)
    P <- matrix( NA , nrow=TP , ncol=ncat)
    P[,1] <- 1
    for (cc in 2:ncat){
        P[,cc] <- exp( Theta \%*\% b )
    }
    P <- P / rowSums(P)    
    return(P)
}
par <- seq( -1.5 , 1.5 , length=K )
names(par) <- paste0("b",1:K)
item_lc <- sirt::xxirt_createDiscItem( name = "LC", par=par, 
                 est=rep(TRUE,K) , P=P_lc )
customItems <- list( item_lc )

# create parameter table
itemtype <- rep( "LC" , 12 )
partable <- sirt::xxirt_createParTable( dat, itemtype=itemtype, customItems=customItems)
partable

#*** estimate model
mod4 <- sirt::xxirt( dat = dat, Theta=Theta, partable = partable, customItems = customItems,
                customTheta= customTheta)
summary(mod4)                
# class probabilities
mod4$probs_Theta
# item response functions
imod4 <- IRT.irfprob( mod5 )
round( imod4[,2,] , 3 )

#****************************************************************************
#******* Model 5: Ordered latent class model

K <- 3 # number of classes
Theta <- diag(K)
Theta <- apply( Theta , 1 , cumsum )

#*** Theta distribution
P_Theta1 <- function( par , Theta , G  ){
    logitprobs <- par[1:(K-1)]
    l1 <- exp( c( logitprobs , 0 ) )    
    probs <- matrix( l1/sum(l1) , ncol=1)    
    return(probs)
}              
par_Theta <- stats::qlogis( rep( 1/K , K-1 ) )
names(par_Theta) <- paste0("pi",1:(K-1) )
customTheta  <- sirt::xxirt_createThetaDistribution( par=par_Theta , 
                est= rep(TRUE,K-1) , P=P_Theta1  )

#*** IRF ordered latent class
P_olc <- function( par, Theta , ncat){
    b <- par
    TP <- nrow(Theta)
    P <- matrix( NA , nrow=TP , ncol=ncat)
    P[,1] <- 1
    for (cc in 2:ncat){
        P[,cc] <- exp( Theta \%*\% b )
    }
    P <- P / rowSums(P)    
    return(P)
}

par <- c( -1 , rep( .5 , , length=K-1 ) )
names(par) <- paste0("b",1:K)
item_olc <- sirt::xxirt_createDiscItem( name = "OLC" , par = par , est = rep(TRUE,K) , 
                    P = P_olc , lower=c( -Inf , 0 , 0 ) )                 
customItems <- list( item_olc )
itemtype <- rep( "OLC" , 12 )
partable <- sirt::xxirt_createParTable( dat, itemtype = itemtype, customItems = customItems)
partable

#*** estimate model
mod5 <- sirt::xxirt( dat=dat, Theta=Theta, partable = partable, customItems = customItems,
                customTheta= customTheta )
summary(mod5)                
# estimated item response functions
imod5 <- IRT.irfprob( mod5 )
round( imod5[,2,] , 3 )

#############################################################################
## EXAMPLE 2: Multiple group models with xxirt
#############################################################################

data(data.math)
dat <- data.math$data
items <- grep( "M[A-Z]" , colnames(dat) , value=TRUE )
I <- length(items)

Theta <- matrix( seq(-8,8,len=31) , ncol=1 )

#****************************************************************************
#******* Model 1: Rasch model, single group

#*** Theta distribution
P_Theta1 <- function( par , Theta , G  ){
    mu <- par[1]
    sigma <- max( par[2] , .01 )
    p1 <- stats::dnorm( Theta[,1] , mean = mu , sd = sigma)
    p1 <- p1 / sum(p1)          
    probs <- matrix( p1 , ncol=1)    
    return(probs)
}
               
par_Theta <- c(0,1)
names(par_Theta) <- c("mu","sigma")
customTheta  <- sirt::xxirt_createThetaDistribution( par=par_Theta , 
                   est= c(FALSE,TRUE) , P=P_Theta1  )
customTheta

#*** IRF 1PL logit
P_1PL <- function( par, Theta , ncat){
    b <- par
    TP <- nrow(Theta)
    P <- matrix( NA , nrow=TP , ncol=ncat)
    P[,2] <- plogis( Theta - b )
    P[,1] <- 1 - P[,2]
    return(P)
}
par <- c("b"=0)
item_1PL <- sirt::xxirt_createDiscItem( name = "1PL", par = par, est = c(TRUE), P = P_1PL)
customItems <- list( item_1PL )

itemtype <- rep( "1PL" , I )
partable <- sirt::xxirt_createParTable( dat[,items] , itemtype = itemtype , 
                       customItems = customItems )
partable <- sirt::xxirt_modifyParTable( partable=partable , parname = "b" , 
                  value = - stats::qlogis( colMeans(dat[,items] )) )

#*** estimate model
mod1 <- sirt::xxirt( dat = dat[,items] , Theta=Theta , partable = partable, 
                customItems = customItems, customTheta= customTheta )
summary(mod1)                

#****************************************************************************
#******* Model 2: Rasch model, multiple groups

#*** Theta distribution
P_Theta2 <- function( par , Theta , G  ){
    mu1 <- par[1]
    mu2 <- par[2]
    sigma1 <- max( par[3] , .01 )
    sigma2 <- max( par[4] , .01 )    
    TP <- nrow(Theta)
    probs <- matrix( NA , nrow=TP , ncol=G)    
    p1 <- stats::dnorm( Theta[,1] , mean = mu1 , sd = sigma1)
    probs[,1] <- p1 / sum(p1)          
    p1 <- stats::dnorm( Theta[,1] , mean = mu2 , sd = sigma2)
    probs[,2] <- p1 / sum(p1)                  
    return(probs)
}              
par_Theta <- c(0,0,1,1)
names(par_Theta) <- c("mu1","mu2","sigma1","sigma2")
customTheta2  <- sirt::xxirt_createThetaDistribution( par=par_Theta , 
                    est= c(FALSE,TRUE,TRUE,TRUE) , P=P_Theta2  )
customTheta2

#*** estimate model
mod2 <- sirt::xxirt( dat=dat[,items], group=dat$female, Theta=Theta, partable=partable, 
           customItems = customItems , customTheta= customTheta2 , maxit=40) 
summary(mod2)
IRT.compareModels(mod1, mod2)

#*** compare results with TAM package
library(TAM)
mod2b <- TAM::tam.mml( resp=dat[,items] , group = dat$female )
summary(mod2b)
IRT.compareModels(mod1, mod2, mod2b)
}	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
