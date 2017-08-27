%% File Name: gom.em.Rd
%% File Version: 1.62
%% File Last Change: 2017-01-18 18:08:38

\name{gom.em}
\alias{gom.em}
\alias{summary.gom}
\alias{logLik.gom}
\alias{anova.gom}
\alias{IRT.irfprob.gom}
\alias{IRT.likelihood.gom}
\alias{IRT.posterior.gom}
\alias{IRT.modelfit.gom}
\alias{summary.IRT.modelfit.gom}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Discrete (Rasch) Grade of Membership Model
}
\description{
This function estimates the grade of membership model (Erosheva, Fienberg
& Joutard, 2007; also called mixed membership model) by the EM algorithm 
assuming a discrete membership score distribution.
}
\usage{
gom.em(dat, K=NULL, problevels=NULL, model="GOM", theta0.k=seq(-5, 5, len=15), 
    xsi0.k=exp(seq(-6, 3, len=15)), max.increment=0.3, numdiff.parm=0.001, 
    maxdevchange=10^(-5), globconv=0.001, maxiter=1000, msteps=4, mstepconv=0.001,
    progress=TRUE)

\method{summary}{gom}(object,...)

\method{anova}{gom}(object,...)

\method{logLik}{gom}(object,...)

\method{IRT.irfprob}{gom}(object,...)

\method{IRT.likelihood}{gom}(object,...)

\method{IRT.posterior}{gom}(object,...)

\method{IRT.modelfit}{gom}(object,...)

\method{summary}{IRT.modelfit.gom}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame with dichotomous responses
}
  \item{K}{
Number of classes (only applies for \code{model="GOM"})
}
  \item{problevels}{
Vector containing probability levels for membership functions 
(only applies for \code{model="GOM"}). If a specific space of probability
levels should be estimated, then a matrix can be supplied (see Example 1,
Model 2a).
}
  \item{model}{
The type of grade of membership model. The default \code{"GOM"} 
is the nonparametric grade of membership model. The probabilities and membership
functions specifications described in Details are called via \code{"GOMRasch"}.
}
  \item{theta0.k}{
Vector of \eqn{\tilde{\theta}_k} grid (applies only for \code{model="GOMRasch"})
}
  \item{xsi0.k}{
Vector of \eqn{\xi_p} grid (applies only for \code{model="GOMRasch"})
}
  \item{max.increment}{
Maximum increment
}
  \item{numdiff.parm}{
Numerical differentiation parameter
}
  \item{maxdevchange}{
Convergence criterion for change in relative deviance
}
  \item{globconv}{
Global convergence criterion for parameter change
}
  \item{maxiter}{
Maximum number of iterations
}
  \item{msteps}{
Number of iterations within a M step
}
  \item{mstepconv}{
Convergence criterion within a M step
}
  \item{progress}{
Display iteration progress? Default is \code{TRUE}.
}
  \item{object}{
Object of class \code{gom}
}
  \item{\dots}{
Further arguments to be passed
}
}
\details{
The item response model of the grade of membership model 
(Erosheva, Fienberg & Junker, 2002; 
Erosheva, Fienberg & Joutard, 2007) with \eqn{K} classes
for dichotomous correct responses \eqn{X_{pi}}
of person \eqn{p} on item \eqn{i} is as follows (\code{model="GOM"})
\deqn{
    P(X_{pi}=1 | g_{p1}, \ldots , g_{pK} ) = \sum_k \lambda_{ik} g_{pk}
\quad , \quad \sum_{k=1}^K g_{pk} = 1
\quad , \quad 0 \leq g_{pk} \leq 1 
                }

In most applications (e.g. Erosheva et al., 2007), the grade of
memebership function \eqn{\{g_{pk}\}} is assumed to follow a Dirichlet
distribution. In our \code{gom.em} implementation             
the membership function is assumed to be discretely represented 
by a grid \eqn{u=(u_1 , \ldots , u_L)} with entries between 0 and 1
(e.g. \code{seq(0,1,length=5)} with \eqn{L=5}). 
The values \eqn{g_{pk}} of the membership function can then
only take values in \eqn{\{ u_1 , \ldots , u_L \}} with the restriction
\eqn{\sum_k g_{pk} \sum_l \bold{1}(g_{pk} = u_l ) = 1}. 
The grid \eqn{u} is specified by using the argument \code{problevels}.
           
The Rasch grade of membership model (\code{model="GOMRasch"}) poses constraints 
on probabilities \eqn{\lambda_{ik}} and membership functions \eqn{g_{pk}}. 
The membership 
function of person \eqn{p} is parametrized by a location parameter \eqn{\theta_p}
and a variability parameter \eqn{\xi_p}. Each class \eqn{k} is represented by
a location parameter \eqn{\tilde{\theta}_k}. The membership function is defined as
\deqn{ g_{pk} \propto 
\exp \left[ - \frac{ (\theta_p - \tilde{\theta}_k)^2 }{2 \xi_p^2 } \right]
}
   
The person parameter \eqn{\theta_p} indicates the usual 'ability', while
\eqn{\xi_p} describes the individual tendency to change between classes
\eqn{1,\ldots,K} and their corresponding locations
\eqn{\tilde{\theta}_1 , \ldots ,\tilde{\theta}_K}.
The extremal class probabilities \eqn{\lambda_{ik}} follow the Rasch model
\deqn{ \lambda_{ik} = invlogit( \tilde{\theta}_k - b_i  ) =
\frac{ \exp( \tilde{\theta}_k - b_i ) }{ 1 + \exp( \tilde{\theta}_k - b_i ) }}

Putting these assumptions together leads to the model equation
\deqn{
    P(X_{pi}=1 | g_{p1}, \ldots , g_{pK} ) =
    P(X_{pi}=1 | \theta_p , \xi_p  ) =
        \sum_k \frac{ \exp( \tilde{\theta}_k - b_i ) }{ 1 + \exp(\tilde{\theta}_k - b_i ) }
        \cdot \exp \left[ - \frac{ (\theta_p - \tilde{\theta}_k)^2 }{2 \xi_p^2 } \right]
                }

In the extreme case of a very small \eqn{\xi_p = \varepsilon > 0} and
\eqn{\theta_p = \theta_0}, the Rasch model is obtained

\deqn{
    P(X_{pi}=1 | \theta_p , \xi_p  ) =
    P(X_{pi}=1 | \theta_0 , \varepsilon  )  =  
        \frac{ \exp( \theta_0 - b_i ) }{ 1 + \exp( \theta_0 - b_i ) }
                }

See Erosheva et al. (2002), Erosheva (2005, 2006) or Galyart (2015)
for a comparison of grade of membership models with latent trait models 
and latent class models.

The grade of membership model is also published under the name
Bernoulli aspect model, see Bingham, Kaban and Fortelius (2009).
}


\value{
A list with following entries:
\item{deviance}{Deviance}
\item{ic}{Information criteria}
\item{item}{Data frame with item parameters}
\item{person}{Data frame with person parameters}
\item{EAP.rel}{EAP reliability (only applies for \code{model="GOMRasch"})}
\item{MAP}{Maximum aposteriori estimate of the membership function 
	}
\item{classdesc}{Descriptives for class membership}	
\item{lambda}{Estimated response probabilities \eqn{\lambda_{ik}}}
\item{se.lambda}{Standard error for stimated response probabilities 
\eqn{\lambda_{ik}}}
\item{mu}{Mean of the distribution of \eqn{(\theta_p , \xi_p)} 
	(only applies for \code{model="GOMRasch"})}
\item{Sigma}{Covariance matrix of \eqn{(\theta_p , \xi_p)} 
	(only applies for \code{model="GOMRasch"})}
\item{b}{Estimated item difficulties (only applies for \code{model="GOMRasch"})}
\item{se.b}{Standard error of estimated difficulties 
(only applies for \code{model="GOMRasch"})}
\item{f.yi.qk}{Individual likelihood}
\item{f.qk.yi}{Individual posterior}
\item{probs}{Array with response probabilities}
\item{n.ik}{Expected counts}
\item{iter}{Number of iterations}
\item{I}{Number of items}
\item{K}{Number of classes}
\item{TP}{Number of discrete integration points for \eqn{(g_{p1},...,g_{pK})}}
\item{theta.k}{Used grid of membership functions}
\item{\dots}{Further values}
}
\references{
Bingham, E., Kaban, A., & Fortelius, M. (2009). 
The aspect Bernoulli model: multiple causes of presences and absences. 
\emph{Pattern Analysis and Applications}, \bold{12(1)}, 55-78.	
	
Erosheva, E. A. (2005). Comparing latent structures of the grade of membership, 
Rasch, and latent class models. \emph{Psychometrika}, \bold{70}, 619-628.

Erosheva, E. A. (2006). \emph{Latent class representation of the grade of membership 
model}. Seattle: University of Washington.

Erosheva, E. A., Fienberg, S. E., & Junker, B. W. (2002). 
Alternative statistical models and representations for large sparse 
multi-dimensional contingency tables. 
\emph{Annales-Faculte Des Sciences Toulouse Mathematiques}, \bold{11},
485-505.

Erosheva, E. A., Fienberg, S. E., & Joutard, C. (2007).
Describing disability through individual-level mixture models 
for multivariate binary data. \emph{Annals of Applied Statistics},
\bold{1}, 502-537.

Galyardt, A. (2015). 
Interpreting mixed membership models: Implications of Erosheva's representation 
theorem. In E. M. Airoldi, D. Blei, E. A. Erosheva, & S. E. Fienberg (Eds.).
\emph{Handbook of Mixed Membership Models} (pp. 39-65). Chapman & Hall.
}
\author{
Alexander Robitzsch
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For joint maximum likelihood estimation of the grade of membership model
see \code{\link{gom.jml}}.

See also the \pkg{mixedMem} package for estimating mixed membership
models by a variational EM algorithm.

The C code of Erosheva et al. (2007) can be downloaded from
\url{http://projecteuclid.org/euclid.aoas/1196438029#supplemental}.

Code from Manrique-Vallier can be downloaded from
\url{http://pages.iu.edu/~dmanriqu/software.html}.

See \url{http://users.ics.aalto.fi/ella/publications/aspect_bernoulli.m}
for a Matlab implementation of the algorithm in Bingham, Kaban and
Fortelius (2009).
}
\examples{
#############################################################################
# EXAMPLE 1: PISA data mathematics
#############################################################################

data(data.pisaMath)
dat <- data.pisaMath$data
dat <- dat[ , grep("M" , colnames(dat)) ]

#***
# Model 1: Discrete GOM with 3 classes and 5 probability levels
problevels <- seq( 0 , 1 , len=5 )
mod1 <- sirt::gom.em( dat , K=3 , problevels ,  model="GOM"  )            
summary(mod1)

\dontrun{
#***
# Model 2: Discrete GOM with 4 classes and 5 probability levels
problevels <- seq( 0 , 1 , len=5 )
mod2 <- sirt::gom.em( dat , K=4 , problevels ,  model="GOM"  )            
summary(mod2)

# model comparison
smod1 <- IRT.modelfit(mod1)
smod2 <- IRT.modelfit(mod2)
IRT.compareModels(smod1,smod2)

#***
# Model 2a: Estimate discrete GOM with 4 classes and restricted space of probability levels
#  the 2nd, 4th and 6th class correspond to "intermediate stages"
problevels <- scan()
 1  0  0  0
.5 .5  0  0
 0  1  0  0
 0 .5 .5  0
 0  0  1  0
 0  0 .5 .5
 0  0  0  1

problevels <- matrix( problevels, ncol=4 , byrow=TRUE)
mod2a <- sirt::gom.em( dat , K=4 , problevels ,  model="GOM" )            
# probability distribution for latent classes
cbind( mod2a$theta.k , mod2a$pi.k )
  ##        [,1] [,2] [,3] [,4]       [,5]
  ##   [1,]  1.0  0.0  0.0  0.0 0.17214630
  ##   [2,]  0.5  0.5  0.0  0.0 0.04965676
  ##   [3,]  0.0  1.0  0.0  0.0 0.09336660
  ##   [4,]  0.0  0.5  0.5  0.0 0.06555719
  ##   [5,]  0.0  0.0  1.0  0.0 0.27523678
  ##   [6,]  0.0  0.0  0.5  0.5 0.08458620
  ##   [7,]  0.0  0.0  0.0  1.0 0.25945016
}

#***
# Model 3: Rasch GOM
mod3 <- sirt::gom.em( dat , model="GOMRasch" , maxiter=20 )            
summary(mod3)

#***
# Model 4: 'Ordinary' Rasch model
mod4 <- sirt::rasch.mml2( dat )
summary(mod4)

\dontrun{
#############################################################################
# EXAMPLE 2: Grade of membership model with 2 classes
#############################################################################

#********* DATASET 1 *************
# define an ordinary 2 latent class model
set.seed(8765)
I <- 10
prob.class1 <- stats::runif( I , 0 , .35 )
prob.class2 <- stats::runif( I , .70 , .95 )
probs <- cbind( prob.class1 , prob.class2 )

# define classes 
N <- 1000
latent.class <- c( rep( 1 , 1/4*N ) , rep( 2,3/4*N ) )

# simulate item responses
dat <- matrix( NA , nrow=N , ncol=I )
for (ii in 1:I){    
    dat[,ii] <- probs[ ii , latent.class ]
    dat[,ii] <- 1 * ( stats::runif(N) < dat[,ii] )
        }
colnames(dat) <- paste0( "I" , 1:I) 

# Model 1: estimate latent class model
mod1 <- sirt::gom.em(dat, K=2, problevels= c(0,1) , model="GOM" )
summary(mod1)
# Model 2: estimate GOM
mod2 <- sirt::gom.em(dat, K=2, problevels= seq(0,1,0.5) , model="GOM" )
summary(mod2)
# estimated distribution
cbind( mod2$theta.k , mod2$pi.k )
  ##       [,1] [,2]        [,3]
  ##  [1,]  1.0  0.0 0.243925644
  ##  [2,]  0.5  0.5 0.006534278
  ##  [3,]  0.0  1.0 0.749540078

#********* DATASET 2 *************
# define a 2-class model with graded membership
set.seed(8765)
I <- 10
prob.class1 <- stats::runif( I , 0 , .35 )
prob.class2 <- stats::runif( I , .70 , .95 )
prob.class3 <- .5*prob.class1+.5*prob.class2  # probabilities for 'fuzzy class'
probs <- cbind( prob.class1 , prob.class2 , prob.class3)
# define classes 
N <- 1000
latent.class <- c( rep(1,round(1/3*N)),rep(2,round(1/2*N)),rep(3,round(1/6*N)))
# simulate item responses
dat <- matrix( NA , nrow=N , ncol=I )
for (ii in 1:I){
    dat[,ii] <- probs[ ii , latent.class ]
    dat[,ii] <- 1 * ( stats::runif(N) < dat[,ii] )
        }
colnames(dat) <- paste0( "I" , 1:I) 

#** Model 1: estimate latent class model
mod1 <- sirt::gom.em(dat, K=2, problevels= c(0,1) , model="GOM" )
summary(mod1)

#** Model 2: estimate GOM
mod2 <- sirt::gom.em(dat, K=2, problevels= seq(0,1,0.5) , model="GOM" )
summary(mod2)
# inspect distribution
cbind( mod2$theta.k , mod2$pi.k )
  ##       [,1] [,2]      [,3]
  ##  [1,]  1.0  0.0 0.3335666
  ##  [2,]  0.5  0.5 0.1810114
  ##  [3,]  0.0  1.0 0.4854220

#***
# Model2m: estimate discrete GOM in mirt
# define latent classes
Theta <- scan( nlines=1)
   1 0   .5 .5    0 1 
Theta <- matrix( Theta , nrow=3 , ncol=2,byrow=TRUE)
# define mirt model
I <- ncol(dat) 
#*** create customized item response function for mirt model
name <- 'gom'
par <- c("a1" = -1 , "a2" = 1 )
est <- c(TRUE, TRUE)
P.gom <- function(par,Theta,ncat){
    # GOM for two extremal classes
    pext1 <- stats::plogis(par[1])
    pext2 <- stats::plogis(par[2]) 
    P1 <- Theta[,1]*pext1 + Theta[,2]*pext2     
    cbind(1-P1, P1)
}
# create item response function
icc_gom <- mirt::createItem(name, par=par, est=est, P=P.gom)
#** define prior for latent class analysis
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
#*** estimate discrete GOM in mirt package
mod2m <- mirt::mirt(dat, 1, rep( "icc_gom",I) , customItems=list("icc_gom"=icc_gom), 
           technical = list( customTheta=Theta , customPriorFun = lca_prior)  )
# correct number of estimated parameters
mod2m@nest <- as.integer(sum(mod.pars$est) + nrow(Theta)-1 ) 
# extract log-likelihood and compute AIC and BIC
mod2m@logLik
( AIC <- -2*mod2m@logLik+2*mod2m@nest )
( BIC <- -2*mod2m@logLik+log(mod2m@Data$N)*mod2m@nest )
# extract coefficients
( cmod2m <- sirt::mirt.wrapper.coef(mod2m) )
# compare estimated distributions
round( cbind( "sirt"  = mod2$pi.k , "mirt" = mod2m@Prior[[1]] ) , 5 )
  ##           sirt    mirt
  ##   [1,] 0.33357 0.33627
  ##   [2,] 0.18101 0.17789
  ##   [3,] 0.48542 0.48584
# compare estimated item parameters
dfr <- data.frame( "sirt" = mod2$item[,4:5] )
dfr$mirt <- apply(cmod2m$coef[ , c("a1" , "a2") ] , 2 , stats::plogis )
round(dfr,4)
  ##      sirt.lam.Cl1 sirt.lam.Cl2 mirt.a1 mirt.a2
  ##   1        0.1157       0.8935  0.1177  0.8934
  ##   2        0.0790       0.8360  0.0804  0.8360
  ##   3        0.0743       0.8165  0.0760  0.8164
  ##   4        0.0398       0.8093  0.0414  0.8094
  ##   5        0.1273       0.7244  0.1289  0.7243
  ##   [...]
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Grade of membership model}
\keyword{Rasch grade of membership model}
\keyword{summary}
