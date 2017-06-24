\name{prob.guttman}
\alias{prob.guttman}
\alias{summary.prob.guttman}
\alias{logLik.prob.guttman}
\alias{anova.prob.guttman}
\alias{IRT.likelihood.prob.guttman}
\alias{IRT.posterior.prob.guttman}
\alias{IRT.irfprob.prob.guttman}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probabilistic Guttman Model
}
\description{
This function estimates the probabilistic Guttman model which
is a special case of an ordered latent trait model (Hanson, 2000;
Proctor, 1970).
}
\usage{
prob.guttman(dat, pid = NULL, guess.equal = FALSE,  slip.equal = FALSE, 
    itemlevel = NULL, conv1 = 0.001, glob.conv = 0.001, mmliter = 500)
    
\method{summary}{prob.guttman}(object,...)

\method{anova}{prob.guttman}(object,...)

\method{logLik}{prob.guttman}(object,...)

\method{IRT.irfprob}{prob.guttman}(object,...)

\method{IRT.likelihood}{prob.guttman}(object,...)

\method{IRT.posterior}{prob.guttman}(object,...)
}



%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
An \eqn{N \times I} data frame of dichotomous item responses
}
  \item{pid}{
Optional vector of person identifiers
}
  \item{guess.equal}{
Should the same guessing parameters for all the items estimated?
}
  \item{slip.equal}{
Should the same slipping parameters for all the items estimated?
}
  \item{itemlevel}{
A vector of item levels of the Guttman scale for each item. If there
are \eqn{K} different item levels, then the Guttman scale possesses
\eqn{K} ordered trait values.
}
  \item{conv1}{
Convergence criterion for item parameters
}
  \item{glob.conv}{
Global convergence criterion for the deviance
}
  \item{mmliter}{
Maximum number of iterations
}
  \item{object}{
Object of class \code{prob.guttman}
}
  \item{\dots}{
Further arguments to be passed
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
An object of class \code{prob.guttman}
\item{person}{Estimated person parameters}
\item{item}{Estimated item parameters}
\item{theta.k}{Ability levels}
\item{trait}{Estimated trait distribution}
\item{ic}{Information criteria}
\item{deviance}{Deviance}
\item{iter}{Number of iterations}
\item{itemdesign}{Specified allocation of items to trait levels}
}
\references{
Hanson, B. (2000). \emph{IRT parameter estimation using the EM algorithm}. 
Technical Report.

Proctor, C. H. (1970). A probabilistic formulation and statistical analysis 
for Guttman scaling. \emph{Psychometrika}, \bold{35}, 73-78.
}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%A summary can be obtained using \code{\link{summary.prob.guttman}}.
%}


\examples{
#############################################################################
# EXAMPLE 1: Dataset Reading
#############################################################################
data(data.read)
dat <- data.read

#***
# Model 1: estimate probabilistic Guttman model
mod1 <- sirt::prob.guttman( dat )
summary(mod1)

#***
# Model 2: probabilistic Guttman model with equal guessing and slipping parameters
mod2 <- sirt::prob.guttman( dat , guess.equal=TRUE , slip.equal=TRUE)
summary(mod2)

#***
# Model 3: Guttman model with three a priori specified item levels
itemlevel <- rep(1,12)
itemlevel[ c(2,5,8,10,12) ] <- 2
itemlevel[ c(3,4,6) ] <- 3
mod3 <- sirt::prob.guttman( dat , itemlevel=itemlevel )
summary(mod3)

\dontrun{
#***
# Model3m: estimate Model 3 in mirt

library(mirt)
# define four ordered latent classes
Theta <- scan(nlines=1)
   0 0 0    1 0 0   1 1 0   1 1 1    
Theta <- matrix( Theta , nrow=4 , ncol=3,byrow=TRUE)

# define mirt model
I <- ncol(dat)  # I = 12
mirtmodel <- mirt::mirt.model("
        # specify factors for each item level
        C1 = 1,7,9,11
        C2 = 2,5,8,10,12 
        C3 = 3,4,6
        ")
# get initial parameter values
mod.pars <- mirt::mirt(dat, model=mirtmodel ,  pars = "values")
# redefine initial parameter values
mod.pars[ mod.pars$name == "d" ,"value" ]  <- -1
mod.pars[ mod.pars$name \%in\% paste0("a",1:3) & mod.pars$est    ,"value" ]  <- 2
mod.pars
# define prior for latent class analysis
lca_prior <- function(Theta,Etable){
  # number of latent Theta classes
  TP <- nrow(Theta)
  # prior in initial iteration
  if ( is.null(Etable) ){ prior <- rep( 1/TP , TP ) }    
  # process Etable (this is correct for datasets without missing data)
  if ( ! is.null(Etable) ){  
    # sum over correct and incorrect expected responses 
    prior <- ( rowSums(Etable[ , seq(1,2*I,2)]) + rowSums(Etable[,seq(2,2*I,2)]) )/I
                 }
  prior <- prior / sum(prior)  
  return(prior)
}
# estimate model in mirt
mod3m <- mirt::mirt(dat, mirtmodel , pars = mod.pars , verbose=TRUE , 
            technical = list( customTheta=Theta , customPriorFun = lca_prior) )
# correct number of estimated parameters
mod3m@nest <- as.integer(sum(mod.pars$est) + nrow(Theta)-1 ) 
# extract log-likelihood and compute AIC and BIC
mod3m@logLik
( AIC <- -2*mod3m@logLik+2*mod3m@nest )
( BIC <- -2*mod3m@logLik+log(mod3m@Data$N)*mod3m@nest )
# compare with information criteria from prob.guttman
mod3$ic
# model fit in mirt
mirt::M2(mod3m)
# extract coefficients
( cmod3m <- sirt::mirt.wrapper.coef(mod3m) )
# compare estimated distributions
round( cbind( "sirt"  = mod3$trait$prob , "mirt" = mod3m@Prior[[1]] ) , 5 )
  ##           sirt    mirt
  ##   [1,] 0.13709 0.13765
  ##   [2,] 0.30266 0.30303
  ##   [3,] 0.15239 0.15085
  ##   [4,] 0.40786 0.40846
# compare estimated item parameters
ipars <- data.frame( "guess.sirt" = mod3$item$guess , 
                     "guess.mirt" = plogis( cmod3m$coef$d ) )
ipars$slip.sirt <- mod3$item$slip
ipars$slip.mirt <- 1-plogis( rowSums(cmod3m$coef[ , c("a1","a2","a3","d") ] ) )
round( ipars , 4 )
  ##      guess.sirt guess.mirt slip.sirt slip.mirt
  ##   1      0.7810     0.7804    0.1383    0.1382
  ##   2      0.4513     0.4517    0.0373    0.0368
  ##   3      0.3203     0.3200    0.0747    0.0751
  ##   4      0.3009     0.3007    0.3082    0.3087
  ##   5      0.5776     0.5779    0.1800    0.1798
  ##   6      0.3758     0.3759    0.3047    0.3051
  ##   7      0.7262     0.7259    0.0625    0.0623
  ##   [...]

#***
# Model 4: Monotone item response function estimated in mirt

# define four ordered latent classes
Theta <- scan(nlines=1)
   0 0 0    1 0 0   1 1 0   1 1 1    
Theta <- matrix( Theta , nrow=4 , ncol=3,byrow=TRUE)

# define mirt model
I <- ncol(dat)  # I = 12
mirtmodel <- mirt::mirt.model("
        # specify factors for each item level
        C1 = 1-12  
        C2 = 1-12
        C3 = 1-12
        ")
# get initial parameter values
mod.pars <- mirt::mirt(dat, model=mirtmodel ,  pars = "values")
# redefine initial parameter values
mod.pars[ mod.pars$name == "d" ,"value" ]  <- -1
mod.pars[ mod.pars$name \%in\% paste0("a",1:3) & mod.pars$est    ,"value" ]  <- .6
# set lower bound to zero ton ensure monotonicity
mod.pars[ mod.pars$name \%in\% paste0("a",1:3)     ,"lbound" ]  <- 0
mod.pars
# estimate model in mirt
mod4 <- mirt::mirt(dat, mirtmodel , pars = mod.pars , verbose=TRUE , 
            technical = list( customTheta=Theta , customPriorFun = lca_prior) )
# correct number of estimated parameters
mod4@nest <- as.integer(sum(mod.pars$est) + nrow(Theta)-1 ) 
# extract coefficients
cmod4 <- sirt::mirt.wrapper.coef(mod4)
cmod4
# compute item response functions
cmod4c <- cmod4$coef[ , c("d" , "a1" , "a2" , "a3" ) ]
probs4 <- t( apply( cmod4c , 1 , FUN = function(ll){ 
                 plogis(cumsum(as.numeric(ll))) } ) )
matplot( 1:4 ,  t(probs4) , type="b" , pch=1:I)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Probabilistic Guttman model}
\keyword{summary}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
