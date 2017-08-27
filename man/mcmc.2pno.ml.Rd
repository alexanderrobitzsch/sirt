%% File Name: mcmc.2pno.ml.Rd
%% File Version: 0.34
%% File Last Change: 2017-01-18 18:08:39

\name{mcmc.2pno.ml}
\alias{mcmc.2pno.ml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Random Item Response Model / Multilevel IRT Model
}
\description{
This function enables the estimation of random item models and multilevel
(or hierarchical) IRT models (Chaimongkol, Huffer & Kamata, 2007;
Fox & Verhagen, 2010; van den Noortgate, de Boeck & Meulders, 2003;
Asparouhov & Muthen, 2012; Muthen & Asparouhov, 2013, 2014). 
Dichotomous response data is supported using a probit link. Normally
distributed responses can also be analyzed. 
See Details for a description of the implemented item
response models.
}
\usage{
mcmc.2pno.ml(dat, group, link="logit" , est.b.M = "h", est.b.Var = "n", 
    est.a.M = "f", est.a.Var = "n", burnin = 500, iter = 1000, 
    N.sampvalues = 1000, progress.iter = 50, prior.sigma2 = c(1, 0.4), 
    prior.sigma.b = c(1, 1), prior.sigma.a = c(1, 1), prior.omega.b = c(1, 1), 
    prior.omega.a = c(1, 0.4) , sigma.b.init=.3 )
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\arguments{
  \item{dat}{
  Data frame with item responses.
}
  \item{group}{
Vector of group identifiers (e.g. classes, schools or countries)
}
\item{link}{Link function. Choices are \code{"logit"} for dichotomous data
	and \code{"normal"} for data under normal distribution assumptions
		}
  \item{est.b.M}{
Estimation type of \eqn{b_i} parameters: \cr \code{n} - non-hierarchical prior
distribution, i.e. \eqn{\omega_b} is set to a very high value and is not 
estimated \cr \code{h} - hierarchical prior distribution with estimated
distribution parameters \eqn{\mu_b} and \eqn{\omega_b}
}
 \item{est.b.Var}{
Estimation type of standard deviations of item difficulties \eqn{b_i}. \cr
\code{n} -- no estimation of the item variance, i.e. \eqn{\sigma_{b,i}} is 
assumed to be zero \cr
\code{i} -- item-specific standard deviation of item difficulties \cr
\code{j} -- a joint standard deviation of all item difficulties is estimated,
i.e. \eqn{\sigma_{b,1} = \ldots  = \sigma_{b,I} = \sigma_b}
}
\item{est.a.M}{
Estimation type of \eqn{a_i} parameters: \cr 
\code{f} - no estimation of item slopes, i.e all item slopes \eqn{a_i} are
fixed at one \cr
\code{n} - non-hierarchical prior distribution, i.e. \eqn{\omega_a =0} \cr 
\code{h} - hierarchical prior distribution with estimated
distribution parameter \eqn{\omega_a}
}
  \item{est.a.Var}{
Estimation type of standard deviations of item slopes \eqn{a_i}. \cr
\code{n} -- no estimation of the item variance \cr
\code{i} -- item-specific standard deviation of item slopes \cr
\code{j} -- a joint standard deviation of all item slopes is estimated,
i.e. \eqn{\sigma_{a,1} = \ldots  = \sigma_{a,I} = \sigma_a}
}
  \item{burnin}{
Number of burnin iterations
}
\item{iter}{
    Total number of iterations
}
  \item{N.sampvalues}{
   Maximum number of sampled values to save
}
\item{progress.iter}{
   Display progress every \code{progress.iter}-th iteration. If no progress
   display is wanted, then choose \code{progress.iter} larger than \code{iter}.
}
  \item{prior.sigma2}{
Prior for Level 2 standard deviation \eqn{\sigma_{L2}}
}
  \item{prior.sigma.b}{
Priors for item difficulty standard deviations \eqn{\sigma_{b,i}}
}
  \item{prior.sigma.a}{
Priors for item difficulty standard deviations \eqn{\sigma_{a,i}}
}
  \item{prior.omega.b}{
Prior for \eqn{\omega_b}
  }
  \item{prior.omega.a}{
Prior for \eqn{\omega_a}
}
\item{sigma.b.init}{
Initial standard deviation for \eqn{\sigma_{b,i}} parameters	
}
}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETAILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\details{
For dichotomous item responses (\code{link="logit"}) of persons \eqn{p} in 
group \eqn{j} on 
item \eqn{i}, the probability of a correct response is defined as
\deqn{P( X_{pji} = 1 | \theta_{pj} ) = \Phi ( a_{ij} \theta_{pj} - b_{ij} )}
The ability \eqn{\theta_{pj}} is decomposed into a Level 1 and a Level 2
effect \deqn{\theta_{pj} = u_j + e_{pj} \quad , \quad 
	u_j \sim N ( 0 , \sigma_{L2}^2 ) \quad , \quad
	e_{pj} \sim N ( 0 , \sigma_{L1}^2 ) }
In a multilevel IRT model (or a random item model), item parameters are
allowed to vary across groups:
\deqn{ b_{ij} \sim N( b_i , \sigma^2_{b,i} ) \quad , \quad
	a_{ij} \sim N( a_i , \sigma^2_{a,i} ) }
In a hierarchical IRT model, a hierarchical distribution of the (main)
item parameters is assumed
\deqn{ b_{i} \sim N( \mu_b , \omega^2_{b} ) \quad , \quad
	a_{i} \sim N( 1 , \omega^2_{a} ) }
Note that for identification purposes, the mean of all item slopes \eqn{a_i}
is set to one. Using the arguments \code{est.b.M}, \code{est.b.Var},
\code{est.a.M} and \code{est.a.Var} defines which variance components 
should be estimated.

For normally distributed item responses (\code{link="normal"}), the model
equations remain the same except the item response model which is now written as
\deqn{ X_{pji} = a_{ij} \theta_{pj} - b_{ij} + \varepsilon_{pji} \quad ,
\quad \varepsilon_{pji} \sim N( 0 , \sigma^2_{res,i} ) }
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\value{
A list of class \code{mcmc.sirt} with following entries:
  \item{mcmcobj}{Object of class \code{mcmc.list}}
  \item{summary.mcmcobj}{Summary of the \code{mcmcobj} object. In this
  summary the Rhat statistic and the mode estimate MAP is included.
  The variable \code{PercSEratio} indicates the proportion of the Monte Carlo
  standard error in relation to the total standard deviation of the
  posterior distribution.}
  \item{ic}{Information criteria (DIC)}
  \item{burnin}{Number of burnin iterations}
  \item{iter}{Total number of iterations}
  \item{theta.chain}{Sampled values of \eqn{\theta_{pj}} parameters}
  \item{theta.chain}{Sampled values of \eqn{u_{j}} parameters}  
  \item{deviance.chain}{Sampled values of Deviance values}  
  \item{EAP.rel}{EAP reliability}
  \item{person}{Data frame with EAP person parameter estimates for 
  	  \eqn{\theta_pj} and their corresponding posterior standard
  	deviations}  
  \item{dat}{Used data frame}
  \item{\dots}{Further values}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\references{
Asparouhov, T. & Muthen, B. (2012). General random effect latent variable 
modeling: Random subjects, items, contexts, and parameters. 
\url{http://www.statmodel.com/papers_date.shtml}.
	
Chaimongkol, S., Huffer, F. W., & Kamata, A. (2007). 
An explanatory differential item functioning (DIF) 
model by the WinBUGS 1.4. 
\emph{Songklanakarin Journal of Science and Technology}, \bold{29}, 449-458.	
	
Fox, J.-P., & Verhagen, A.-J. (2010). Random item effects modeling for 
cross-national survey data. 
In E. Davidov, P. Schmidt, & J. Billiet (Eds.), 
\emph{Cross-cultural Analysis: Methods and Applications} 
(pp. 467-488), London: Routledge Academic.

Muthen, B. & Asparouhov, T. (2013). New methods for the study of measurement 
invariance with many groups. \url{http://www.statmodel.com/papers_date.shtml}

Muthen, B. & Asparouhov, T. (2014). 
Item response modeling in Mplus: A multi-dimensional, multi-level, 
and multi-timepoint example. In W. Linden & R. Hambleton (2014). 
\emph{Handbook of item response theory: Models, statistical tools, and 
applications}. \url{http://www.statmodel.com/papers_date.shtml}

van den Noortgate, W., De Boeck, P., & Meulders, M. (2003). 
Cross-classification multilevel logistic models in psychometrics. 
\emph{Journal of Educational and Behavioral Statistics}, \bold{28}, 369-386.
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
	
For MCMC estimation of three-parameter (testlet) models see 
\code{\link{mcmc.3pno.testlet}}.	

See also the \pkg{MLIRT} package (http://www.jean-paulfox.com).

For more flexible estimation of multilevel IRT models see the
\pkg{MCMCglmm} and \pkg{lme4} packages.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\examples{
\dontrun{	
#############################################################################
# EXAMPLE 1: Dataset Multilevel data.ml1 - dichotomous items
#############################################################################
data(data.ml1)
dat <- data.ml1[,-1]
group <- data.ml1$group
# just for a try use a very small number of iterations
burnin <- 50 ; iter <- 100

#***
# Model 1: 1PNO with no cluster item effects
mod1 <- sirt::mcmc.2pno.ml( dat , group , est.b.Var="n" , burnin=burnin , iter=iter )
summary(mod1)	# summary
plot(mod1,layout=2,ask=TRUE) # plot results
# write results to coda file
mcmclist2coda( mod1$mcmcobj , name = "data.ml1_mod1" )

#***
# Model 2: 1PNO with cluster item effects of item difficulties 
mod2 <- sirt::mcmc.2pno.ml( dat , group , est.b.Var="i" , burnin=burnin , iter=iter )
summary(mod2)
plot(mod2, ask=TRUE , layout=2 )

#***
# Model 3: 2PNO with cluster item effects of item difficulties but
#          joint item slopes
mod3 <- sirt::mcmc.2pno.ml( dat , group , est.b.Var="i" , est.a.M="h" , 
              burnin=burnin , iter=iter )
summary(mod3)

#***
# Model 4: 2PNO with cluster item effects of item difficulties and
#          cluster item effects with a jointly estimated SD
mod4 <- sirt::mcmc.2pno.ml( dat , group , est.b.Var="i" , est.a.M="h" ,
                est.a.Var="j" , burnin=burnin , iter=iter )
summary(mod4)

#############################################################################
# EXAMPLE 2: Dataset Multilevel data.ml2 - polytomous items
#            assuming a normal distribution for polytomous items
#############################################################################
data(data.ml2)
dat <- data.ml2[,-1]
group <- data.ml2$group
# set iterations for all examples (too few!!)
burnin <- 100 ; iter <- 500

#***
# Model 1: no intercept variance, no slopes
mod1 <- sirt::mcmc.2pno.ml( dat=dat , group=group , est.b.Var="n" , 
             burnin=burnin , iter=iter , link="normal" ,  progress.iter=20  )
summary(mod1)

#***
# Model 2a: itemwise intercept variance, no slopes
mod2a <- sirt::mcmc.2pno.ml( dat=dat , group=group , est.b.Var="i" , 
            burnin=burnin , iter=iter ,link="normal" ,  progress.iter=20  )
summary(mod2a)

#***
# Model 2b: homogeneous intercept variance, no slopes
mod2b <- sirt::mcmc.2pno.ml( dat=dat , group=group , est.b.Var="j" , 
              burnin=burnin , iter=iter ,link="normal" ,  progress.iter=20  )
summary(mod2b)

#***
# Model 3: intercept variance and slope variances
#          hierarchical item and slope parameters
mod3 <- sirt::mcmc.2pno.ml( dat=dat , group=group , 
               est.b.M="h" , est.b.Var="i" , est.a.M="h" , est.a.Var="i" , 
               burnin=burnin , iter=iter ,link="normal" ,  progress.iter=20  )
summary(mod3)

#############################################################################
# EXAMPLE 3: Simulated random effects model | dichotomous items
#############################################################################
set.seed(7698)

#*** model parameters
sig2.lev2 <- .3^2   # theta level 2 variance
sig2.lev1 <- .8^2   # theta level 1 variance 
G <- 100            # number of groups
n <- 20             # number of persons within a group
I <- 12             # number of items
#*** simuate theta
theta2 <- stats::rnorm( G , sd = sqrt(sig2.lev2) )
theta1 <- stats::rnorm( n*G , sd =  sqrt(sig2.lev1) )
theta  <- theta1 + rep( theta2 , each=n )
#*** item difficulties
b <- seq( -2 , 2 , len=I )
#*** define group identifier
group <- 1000 + rep(1:G , each=n )
#*** SD of group specific difficulties for items 3 and 5
sigma.item <- rep(0,I)
sigma.item[c(3,5)] <- 1
#*** simulate group specific item difficulties
b.class <- sapply( sigma.item , FUN = function(sii){ stats::rnorm( G , sd = sii ) } )
b.class <- b.class[ rep( 1:G ,each=n ) , ]
b <- matrix( b , n*G , I , byrow=TRUE ) + b.class
#*** simulate item responses
m1 <- stats::pnorm( theta - b ) 
dat <- 1 * ( m1 > matrix( stats::runif( n*G*I ) , n*G , I ) )

#*** estimate model
mod <- sirt::mcmc.2pno.ml( dat , group=group , burnin=burnin , iter=iter , 
            est.b.M="n" , est.b.Var="i" , progress.iter=20)
summary(mod)
plot(mod , layout=2 , ask=TRUE )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Markov Chain Monte Carlo (MCMC)}
\keyword{Multilevel item response model}
\keyword{Multilevel DIF}
