%% File Name: mcmc.3pno.testlet.Rd
%% File Version: 1.17

\name{mcmc.3pno.testlet}
\alias{mcmc.3pno.testlet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
3PNO Testlet Model
}
\description{
This function estimates the 3PNO testlet model (Wang, Bradlow & Wainer, 2002, 2007)
by Markov Chain Monte Carlo methods (Glas, 2012).
}
\usage{
mcmc.3pno.testlet(dat, testlets = rep(NA, ncol(dat)), 
   weights = NULL, est.slope = TRUE, est.guess = TRUE, guess.prior = NULL, 
   testlet.variance.prior = c(1, 0.2), burnin = 500, iter = 1000, 
   N.sampvalues = 1000, progress.iter = 50, save.theta = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame with dichotomous item responses for \eqn{N} persons and \eqn{I} items
}
  \item{testlets}{
An integer or character vector which indicates the allocation of items to 
testlets. Same entries corresponds to same testlets. 
If an entry is \code{NA}, then this item does not belong to any testlet. 
}
  \item{weights}{
An optional vector with student sample weights
}
  \item{est.slope}{
Should item slopes be estimated? The default is \code{TRUE}.
}
  \item{est.guess}{
Should guessing parameters be estimated? The default is \code{TRUE}.
}
  \item{guess.prior}{
A vector of length two or a matrix with \eqn{I} items and two columns 
which defines the beta prior distribution of guessing
parameters. The default is a non-informative prior, i.e. the Beta(1,1)
distribution.
}
  \item{testlet.variance.prior}{
A vector of length two which defines the (joint) prior for testlet variances
assuming an inverse chi-squared distribution.
The first entry is the effective sample size of the prior while the second
entry defines the prior variance of the testlet. The default of \code{c(1,.2)}
means that the prior sample size is 1 and the prior testlet variance is .2.
}
  \item{burnin}{
Number of burnin iterations
}
  \item{iter}{
Number of iterations
}
  \item{N.sampvalues}{
Maximum number of sampled values to save
}
  \item{progress.iter}{
Display progress every \code{progress.iter}-th iteration. If no progress
display is wanted, then choose \code{progress.iter} larger than \code{iter}.
}
  \item{save.theta}{
Should theta values be saved?
}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETAILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\details{
The testlet response model for person \eqn{p} at item \eqn{i} 
is defined as 
\deqn{ P(X_{pi} = 1 ) = c_i + ( 1 - c_i ) 
        \Phi ( a_i \theta_p + \gamma_{p,t(i)} + b_i    ) \quad , \quad
    \theta_p \sim N ( 0  ,1 ) , \gamma_{p,t(i)} \sim N( 0 , \sigma^2_t ) }

In case of \code{est.slope=FALSE}, all item slopes \eqn{a_i} are set to 1. Then
a variance \eqn{\sigma^2} of the \eqn{\theta_p} distribution is estimated
which is called the Rasch testlet model in the literature (Wang & Wilson, 2005).

In case of \code{est.guess=FALSE}, all guessing parameters \eqn{c_i} are
set to 0.

After fitting the testlet model, marginal item parameters are calculated (integrating
out testlet effects \eqn{\gamma_{p,t(i)}}) according the defining response equation
\deqn{ P(X_{pi} = 1 ) = c_i + ( 1 - c_i ) 
        \Phi ( a_i^\ast \theta_p + b_i^\ast    ) }

}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\value{
A list of class \code{mcmc.sirt} with following entries:
  \item{mcmcobj}{Object of class \code{mcmc.list} containing item parameters
  (\code{b_marg} and \code{a_marg} denote marginal item parameters)
  and person parameters (if requested)}
  \item{summary.mcmcobj}{Summary of the \code{mcmcobj} object. In this
  summary the Rhat statistic and the mode estimate MAP is included.
  The variable \code{PercSEratio} indicates the proportion of the Monte Carlo
  standard error in relation to the total standard deviation of the
  posterior distribution.}
  \item{ic}{Information criteria (DIC)}
  \item{burnin}{Number of burnin iterations}
  \item{iter}{Total number of iterations}
  \item{theta.chain}{Sampled values of \eqn{\theta_p} parameters}
  \item{deviance.chain}{Sampled values of deviance values}  
  \item{EAP.rel}{EAP reliability}
  \item{person}{Data frame with EAP person parameter estimates for 
        \eqn{\theta_p} and their corresponding posterior standard
      deviations and for all testlet effects}  
  \item{dat}{Used data frame}
  \item{weights}{Used student weights}
  \item{\dots}{Further values}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\references{
Glas, C. A. W. (2012). \emph{Estimating and testing the extended testlet model.}
LSAC Research Report Series, RR 12-03.

Wainer, H., Bradlow, E. T., & Wang, X. (2007). 
\emph{Testlet response theory and its applications}. 
Cambridge: Cambridge University Press.

Wang, W.-C., & Wilson, M. (2005). The Rasch testlet model.
\emph{Applied Psychological Measurement}, \bold{29}, 126-149. 

Wang, X., Bradlow, E. T., & Wainer, H. (2002). A general Bayesian model 
for testlets: Theory and applications. 
\emph{Applied Psychological Measurement}, \bold{26}, 109-128.
}

\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
S3 methods: \code{\link{summary.mcmc.sirt}}, \code{\link{plot.mcmc.sirt}}    
}

%For estimating testlet models using the \pkg{lme4} package see
%\code{\link{rasch.testlet.glmer}}.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\examples{
\dontrun{    
#############################################################################
# EXAMPLE 1: Dataset Reading
#############################################################################
data(data.read)
dat <- data.read
I <- ncol(dat)

# set burnin and total number of iterations here (CHANGE THIS!)
burnin <- 200
iter <- 500

#***
# Model 1: 1PNO model
mod1 <- sirt::mcmc.3pno.testlet( dat ,  est.slope=FALSE , est.guess=FALSE , 
            burnin=burnin, iter=iter )
summary(mod1)
plot(mod1,ask=TRUE) # plot MCMC chains in coda style
plot(mod1,ask=TRUE , layout=2) # plot MCMC output in different layout

#***
# Model 2: 3PNO model with Beta(5,17) prior for guessing parameters
mod2 <- sirt::mcmc.3pno.testlet( dat ,  guess.prior=c(5,17) ,
               burnin=burnin, iter=iter )
summary(mod2)

#***
# Model 3: Rasch (1PNO) testlet model
testlets <- substring( colnames(dat) , 1 , 1 )
mod3 <- sirt::mcmc.3pno.testlet( dat ,  testlets=testlets ,  est.slope=FALSE , 
           est.guess=FALSE , burnin=burnin, iter=iter )
summary(mod3)

#***
# Model 4: 3PNO testlet model with (almost) fixed guessing parameters .25
mod4 <- sirt::mcmc.3pno.testlet( dat ,  guess.prior=1000*c(25,75) , testlets=testlets , 
              burnin=burnin, iter=iter )
summary(mod4)
plot(mod4, ask=TRUE, layout=2)

#############################################################################
# EXAMPLE 2: Simulated data according to the Rasch testlet model
#############################################################################
set.seed(678)

N <- 3000   # number of persons
I <- 4      # number of items per testlet
TT <- 3     # number of testlets

ITT <- I*TT
b <- round( stats::rnorm( ITT , mean=0 , sd = 1 ) , 2 )
sd0 <- 1 # sd trait
sdt <- seq( 0 , 2 , len=TT ) # sd testlets

# simulate theta
theta <- stats::rnorm( N , sd = sd0 )
# simulate testlets
ut <- matrix(0,nrow=N , ncol=TT )
for (tt in 1:TT){ 
    ut[,tt] <- stats::rnorm( N , sd = sdt[tt] ) 
}
ut <- ut[ , rep(1:TT,each=I) ]
# calculate response probability
prob <- matrix( stats::pnorm( theta + ut + matrix( b , nrow=N , ncol=ITT , 
            byrow=TRUE ) ) , N, ITT)
Y <- (matrix( stats::runif(N*ITT) , N , ITT) < prob )*1
colMeans(Y)

# define testlets
testlets <- rep(1:TT , each=I )

burnin <- 300
iter <- 1000

#***
# Model 1: 1PNO model (without testlet structure)
mod1 <- sirt::mcmc.3pno.testlet( dat=Y ,  est.slope=FALSE , est.guess=FALSE , 
            burnin=burnin, iter=iter , testlets= testlets )        
summary(mod1)        

summ1 <- mod1$summary.mcmcobj
# compare item parameters
cbind( b , summ1[ grep("b" , summ1$parameter ) , "Mean" ] )
# Testlet standard deviations
cbind( sdt , summ1[ grep("sigma\\.testlet" , summ1$parameter ) , "Mean" ] )

#***
# Model 2: 1PNO model (without testlet structure)
mod2 <- sirt::mcmc.3pno.testlet( dat=Y ,  est.slope=TRUE , est.guess=FALSE , 
           burnin=burnin, iter=iter , testlets= testlets )        
summary(mod2)        

summ2 <- mod2$summary.mcmcobj
# compare item parameters
cbind( b , summ2[ grep("b\\[" , summ2$parameter ) , "Mean" ] )
# item discriminations
cbind( sd0 , summ2[ grep("a\\[" , summ2$parameter ) , "Mean" ] )
# Testlet standard deviations
cbind( sdt , summ2[ grep("sigma\\.testlet" , summ2$parameter ) , "Mean" ] )

#############################################################################
# EXAMPLE 3: Simulated data according to the 2PNO testlet model
#############################################################################
set.seed(678)

N <- 3000    # number of persons
I <- 3      # number of items per testlet
TT <- 5    # number of testlets

ITT <- I*TT
b <- round( stats::rnorm( ITT , mean=0 , sd = 1 ) , 2 )
a <- round( stats::runif( ITT , 0.5 , 2 ) ,2)
sdt <- seq( 0 , 2 , len=TT ) # sd testlets
sd0 <- 1

# simulate theta
theta <- stats::rnorm( N , sd = sd0 )
# simulate testlets
ut <- matrix(0,nrow=N , ncol=TT )
for (tt in 1:TT){ 
   ut[,tt] <- stats::rnorm( N , sd = sdt[tt] ) 
}
ut <- ut[ , rep(1:TT,each=I) ]
# calculate response probability
bM <- matrix( b , nrow=N , ncol=ITT , byrow=TRUE )
aM <- matrix( a , nrow=N , ncol=ITT , byrow=TRUE )
prob <- matrix( stats::pnorm( aM*theta + ut + bM ) , N, ITT)
Y <- (matrix( stats::runif(N*ITT) , N , ITT) < prob )*1
colMeans(Y)

# define testlets
testlets <- rep(1:TT , each=I )

burnin <- 500
iter <- 1500

#***
# Model 1: 2PNO model
mod1 <- sirt::mcmc.3pno.testlet( dat=Y ,  est.slope=TRUE , est.guess=FALSE , 
             burnin=burnin, iter=iter , testlets= testlets )                               
summary(mod1)        

summ1 <- mod1$summary.mcmcobj
# compare item parameters
cbind( b , summ1[ grep("b\\[" , summ1$parameter ) , "Mean" ] )
# item discriminations
cbind( a , summ1[ grep("a\\[" , summ1$parameter ) , "Mean" ] )
# Testlet standard deviations
cbind( sdt , summ1[ grep("sigma\\.testlet" , summ1$parameter ) , "Mean" ] )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Testlet model}
\keyword{Testlets}
\keyword{Markov Chain Monte Carlo (MCMC)}
% __ONLY ONE__ keyword per line
