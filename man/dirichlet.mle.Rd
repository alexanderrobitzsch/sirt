%% File Name: dirichlet.mle.Rd
%% File Version: 0.16

\name{dirichlet.mle}
\alias{dirichlet.mle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum Likelihood Estimation of the Dirichlet Distribution
}
\description{
Maximum likelihood estimation of the parameters of the Dirichlet
distribution
}
\usage{
dirichlet.mle(x, weights=NULL , eps = 10^(-5), convcrit = 1e-05 , maxit=1000, 
     oldfac = .3 , progress=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Data frame with \eqn{N} observations and \eqn{K} variables
of a Dirichlet distribution
}
\item{weights}{Optional vector of frequency weights}
  \item{eps}{
Tolerance number which is added to prevent from logarithms of zero
}
  \item{convcrit}{
Convergence criterion
}
  \item{maxit}{
Maximum number of iterations
}
\item{oldfac}{Convergence acceleration factor. It must be a parameter
between 0 and 1.}
  \item{progress}{
Display iteration progress?
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with following entries
  \item{alpha}{Vector of \eqn{\alpha} parameters}
  \item{alpha0}{The concentration parameter \eqn{\alpha_0 = \sum_k \alpha_k} }
  \item{xsi}{Vector of proportions \eqn{\xi_k = \alpha_k / \alpha_0} } 
}
\references{
Minka, T. P. (2012). \emph{Estimating a Dirichlet distribution}.
Technical Report.
}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For simulating Dirichlet vectors with matrix-wise
\eqn{\bold{\alpha}} parameters see \code{\link{dirichlet.simul}}.

For a variety of functions concerning the Dirichlet distribution
see the \pkg{DirichletReg} package.
}

%% ~~objects to See Also as \code{\link{help}}, ~~~

\examples{
#############################################################################
# EXAMPLE 1: Simulate and estimate Dirichlet distribution
#############################################################################	
	
# (1) simulate data
set.seed(789)
N <- 200
probs <- c(.5 , .3 , .2 )
alpha0 <- .5
alpha <- alpha0*probs
alpha <- matrix( alpha , nrow=N , ncol=length(alpha) , byrow=TRUE  )
x <- sirt::dirichlet.simul( alpha )

# (2) estimate Dirichlet parameters
dirichlet.mle(x)
  ##   $alpha
  ##   [1] 0.24507708 0.14470944 0.09590745
  ##   $alpha0
  ##   [1] 0.485694
  ##   $xsi
  ##   [1] 0.5045916 0.2979437 0.1974648

\dontrun{  
#############################################################################
# EXAMPLE 2: Fitting Dirichlet distribution with frequency weights
#############################################################################

# define observed data
x <- scan( nlines=1)
    1 0   0 1   .5 .5
x <- matrix( x , nrow=3 , ncol=2 , byrow=TRUE)    

# transform observations x into (0,1)
eps <- .01
x <- ( x + eps ) / ( 1 + 2 * eps )

# compare results with likelihood fitting package maxLik
miceadds::library_install("maxLik")
# define likelihood function
dirichlet.ll <- function(param) {
    ll <- sum( weights * log( ddirichlet( x , param ) ) )
    ll
}

#*** weights 10-10-1
weights <- c(10, 10 , 1 )
mod1a <- sirt::dirichlet.mle( x , weights= weights )
mod1a
# estimation in maxLik
mod1b <- maxLik::maxLik(loglik, start=c(.5,.5)) 
print( mod1b )
coef( mod1b )

#*** weights 10-10-10
weights <- c(10, 10 , 10 )
mod2a <- sirt::dirichlet.mle( x , weights= weights )
mod2a
# estimation in maxLik
mod2b <- maxLik::maxLik(loglik, start=c(.5,.5)) 
print( mod2b )
coef( mod2b )

#*** weights 30-10-2
weights <- c(30, 10 , 2 )
mod3a <- sirt::dirichlet.mle( x , weights= weights )
mod3a
# estimation in maxLik
mod3b <- maxLik::maxLik(loglik, start=c(.25,.25)) 
print( mod3b )
coef( mod3b )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Dirichlet distribution}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
