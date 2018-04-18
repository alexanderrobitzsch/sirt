%% File Name: sim.rasch.dep.Rd
%% File Version: 0.24

\name{sim.rasch.dep}
\alias{sim.rasch.dep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of the Rasch Model with Locally Dependent Responses
}
\description{
This function simulates dichotomous item responses where for some 
itemclusters residual correlations can be defined.
}
\usage{
sim.rasch.dep(theta, b, itemcluster, rho)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
Vector of person abilities of length \eqn{N}
}
  \item{b}{
Vector of item difficulties of length \eqn{I}
}
  \item{itemcluster}{
Vector of integers (including 0) of length \eqn{I}. Different integers 
correspond to different itemclusters.
}
  \item{rho}{
Vector of residual correlations. The length of vector must be equal
to the number of itemclusters.
}
}


%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
An \eqn{N \times I} data frame of dichotomous item responses.
}
\author{
Alexander Robitzsch
}
\note{
The specification of the simulation models follows a marginal interpretation
of the latent trait. Local dependencies are only interpreted as nuisance
and not of substantive interest. If local dependencies should be substantively 
interpreted, a testlet model seems preferable 
(see \code{\link{mcmc.3pno.testlet}}).
}




%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
To simulate the generalized logistic item response model see 
\code{\link{sim.raschtype}}. Ramsay's quotient model can be simulated 
using \code{\link{sim.qm.ramsay}}.

Marginal item reponse models for locally dependent item responses can be 
estimated with \code{\link{rasch.copula2}}, \code{\link{rasch.pairwise}} or
\code{\link{rasch.pairwise.itemcluster}}.
}
\examples{
#############################################################################
# EXAMPLE 1: 11 Items: 2 itemclusters with 2 resp. 3 dependent items
#             and 6 independent items
#############################################################################	

set.seed(7654)
I <- 11                             # number of items
n <- 1500                           # number of persons
b <- seq(-2,2, len=I)               # item difficulties
theta <- stats::rnorm( n , sd = 1 )        # person abilities
# itemcluster
itemcluster <- rep(0,I)
itemcluster[ c(3,5)] <- 1
itemcluster[c(2,4,9)] <- 2
# residual correlations
rho <- c( .7 , .5 )

# simulate data
dat <- sirt::sim.rasch.dep( theta , b , itemcluster , rho )
colnames(dat) <- paste("I" , seq(1,ncol(dat)) , sep="")

# estimate Rasch copula model
mod1 <- sirt::rasch.copula2( dat , itemcluster = itemcluster )
summary(mod1)

# compare result with Rasch model estimation in rasch.copula
# delta must be set to zero
mod2 <- sirt::rasch.copula2( dat , itemcluster = itemcluster , delta = c(0,0) ,
            est.delta = c(0,0)  )
summary(mod2)

# estimate Rasch model with rasch.mml2 function
mod3 <- sirt::rasch.mml2( dat )
summary(mod3)

\dontrun{
#############################################################################
# EXAMPLE 2: 12 Items: Cluster 1 -> Items 1,...,4; 
#       Cluster 2 -> Items 6,...,9; Cluster 3 -> Items 10,11,12
#############################################################################

set.seed(7896)
I <- 12                             # number of items
n <- 450                            # number of persons
b <- seq(-2,2, len=I)               # item difficulties
b <- sample(b)                      # sample item difficulties
theta <- stats::rnorm( n , sd = 1 )        # person abilities
# itemcluster
itemcluster <- rep(0,I)
itemcluster[ 1:4 ] <- 1
itemcluster[ 6:9 ] <- 2
itemcluster[ 10:12 ] <- 3
# residual correlations
rho <- c( .55 , .25 , .45 )

# simulate data
dat <- sirt::sim.rasch.dep( theta , b , itemcluster , rho )
colnames(dat) <- paste("I" , seq(1,ncol(dat)) , sep="")

# estimate Rasch copula model 
mod1 <- sirt::rasch.copula2( dat , itemcluster = itemcluster , numdiff.parm = .001 )
summary(mod1)

# Rasch model estimation
mod2 <- sirt::rasch.copula2( dat , itemcluster = itemcluster , 
            delta = rep(0,3) , est.delta = rep(0,3) )
summary(mod2)

# estimation with pairwise Rasch model
mod3 <- sirt::rasch.pairwise( dat )
summary(mod3) 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Simulating IRT models}
%\keyword{ ~kwd1 }

