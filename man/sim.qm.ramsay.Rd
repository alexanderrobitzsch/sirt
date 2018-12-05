%% File Name: sim.qm.ramsay.Rd
%% File Version: 0.37

\name{sim.qm.ramsay}
\alias{sim.qm.ramsay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate from Ramsay's Quotient Model
}
\description{
This function simulates dichotomous item response data
according to Ramsay's quotient model (Ramsay, 1989).
}
\usage{
sim.qm.ramsay(theta, b, K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
Vector of of length \eqn{N} person parameters (must be positive!)
}
  \item{b}{
Vector of length \eqn{I} of item difficulties (must be positive)
}
  \item{K}{
Vector of length \eqn{I} of guessing parameters (must be positive)
}
}
\details{
Ramsay's quotient model (Ramsay, 1989) is defined by the equation
\deqn{P(X_{pi}=1 | \theta_p )=\frac{ \exp { ( \theta_p / b_i ) } }
    { K_i + \exp { ( \theta_p / b_i ) } }}
}
\value{
An \eqn{N \times I}  data frame with dichotomous item responses.
}

\references{
Ramsay, J. O. (1989). A comparison of three simple test theory models.
\emph{Psychometrika, 54}, 487-499.

van der Maas, H. J. L., Molenaar, D., Maris, G., Kievit, R. A., &
Borsboom, D. (2011).
Cognitive psychology meets psychometric theory: On the relation between
process models for decision making and latent variable models for
individual differences.
\emph{Psychological Review, 318}, 339-356.
}


%\author{
%Alexander Robitzsch
%}




%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{rasch.mml2}} for estimating Ramsay's quotient model.

See \code{\link{sim.raschtype}} for simulating response data from
the generalized logistic item response model.
}

%% ~~objects to See Also as \code{\link{help}}, ~~~

\examples{
#############################################################################
# EXAMPLE 1: Estimate Ramsay Quotient Model with rasch.mml2
#############################################################################

set.seed(657)
# simulate data according to the Ramsay model
N <- 1000       # persons
I <- 11         # items
theta <- exp( stats::rnorm( N ) )  # person ability
b <- exp( seq(-2,2,len=I))  # item difficulty
K <- rep( 3, I )           # K parameter (=> guessing)

# apply simulation function
dat <- sirt::sim.qm.ramsay( theta, b, K )

#***
# analysis
mmliter <- 50       # maximum number of iterations
I <- ncol(dat)
fixed.K <- rep( 3, I )

# Ramsay QM with fixed K parameter (K=3 in fixed.K specification)
mod1 <- sirt::rasch.mml2( dat, mmliter=mmliter, irtmodel="ramsay.qm",
              fixed.K=fixed.K )
summary(mod1)

# Ramsay QM with joint estimated K parameters
mod2 <- sirt::rasch.mml2( dat, mmliter=mmliter, irtmodel="ramsay.qm",
             est.K=rep(1,I)  )
summary(mod2)

\dontrun{
# Ramsay QM with itemwise estimated K parameters
mod3 <- sirt::rasch.mml2( dat, mmliter=mmliter, irtmodel="ramsay.qm",
              est.K=1:I  )
summary(mod3)

# Rasch model
mod4 <- sirt::rasch.mml2( dat )
summary(mod4)

# generalized logistic model
mod5 <- sirt::rasch.mml2( dat, est.alpha=TRUE, mmliter=mmliter)
summary(mod5)

# 2PL model
mod6 <- sirt::rasch.mml2( dat, est.a=rep(1,I) )
summary(mod6)

# Difficulty + Guessing (b+c) Model
mod7 <- sirt::rasch.mml2( dat, est.c=rep(1,I) )
summary(mod7)

# estimate separate guessing (c) parameters
mod8 <- sirt::rasch.mml2( dat, est.c=1:I  )
summary(mod8)

#*** estimate Model 1 with user defined function in mirt package

# create user defined function for Ramsay's quotient model
name <- 'ramsayqm'
par <- c("K"=3, "b"=1 )
est <- c(TRUE, TRUE)
P.ramsay <- function(par,Theta){
     eps <- .01
     K <- par[1]
     b <- par[2]
     num <- exp( exp( Theta[,1] ) / b )
     denom <- K + num
     P1 <- num / denom
     P1 <- eps + ( 1 - 2*eps ) * P1
     cbind(1-P1, P1)
}

# create item response function
ramsayqm <- mirt::createItem(name, par=par, est=est, P=P.ramsay)
# define parameters to be estimated
mod1m.pars <- mirt::mirt(dat, 1, rep( "ramsayqm",I),
                   customItems=list("ramsayqm"=ramsayqm), pars="values")
mod1m.pars[ mod1m.pars$name=="K", "est" ] <- FALSE
# define Theta design matrix
Theta <- matrix( seq(-3,3,len=10), ncol=1)
# estimate model
mod1m <- mirt::mirt(dat, 1, rep( "ramsayqm",I), customItems=list("ramsayqm"=ramsayqm),
               pars=mod1m.pars, verbose=TRUE,
               technical=list( customTheta=Theta, NCYCLES=50)
                )
print(mod1m)
summary(mod1m)
cmod1m <- sirt::mirt.wrapper.coef( mod1m )$coef
# compare simulated and estimated values
dfr <- cbind( b, cmod1m$b, exp(mod1$item$b ) )
colnames(dfr) <- c("simulated", "mirt", "sirt_rasch.mml2")
round( dfr, 2 )
  ##      simulated mirt sirt_rasch.mml2
  ## [1,]      0.14 0.11            0.11
  ## [2,]      0.20 0.17            0.18
  ## [3,]      0.30 0.27            0.29
  ## [4,]      0.45 0.42            0.43
  ## [5,]      0.67 0.65            0.67
  ## [6,]      1.00 1.00            1.01
  ## [7,]      1.49 1.53            1.54
  ## [8,]      2.23 2.21            2.21
  ## [9,]      3.32 3.00            2.98
  ##[10,]      4.95 5.22            5.09
  ##[11,]      7.39 5.62            5.51
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Simulating IRT models}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
