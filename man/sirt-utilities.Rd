%% File Name: sirt-utilities.Rd
%% File Version: 0.13

\name{sirt-utilities}
\alias{sirt-utilities}
\alias{tracemat}
\alias{pow}
\alias{soft_thresholding}
\alias{hard_thresholding}
\alias{bounds_parameters}
\alias{dimproper}
\alias{ginverse_sym}
\alias{sirt_rbind_fill}
\alias{sirt_fisherz}
\alias{sirt_antifisherz}

\title{Utility Functions in \pkg{sirt}}

\description{
Utility functions in \pkg{sirt}. 
}

\usage{
# bounds entries in a vector
bounds_parameters( pars , lower = NULL , upper = NULL)
	
# improper density function which always returns a value of 1
dimproper(x)

# generalized inverse of a symmetric function
ginverse_sym(A, eps= 1E-8)

# hard thresholding function
hard_thresholding(x, lambda)	

# soft thresholding function
soft_thresholding(x, lambda)
	
# power function x^a, like in Cpp
pow(x, a)

# trace of a matrix
tracemat(A)

# plyr::rbind.fill implementation in sirt
sirt_rbind_fill(x, y)

# Fisher-z transformation, see psych::fisherz
sirt_fisherz(rho)
# inverse Fisher-z transformation, see psych::fisherz2r
sirt_antifisherz(z)
}

\arguments{
\item{pars}{Numeric vector}
\item{lower}{Numeric vector}
\item{upper}{Numeric vector}
	\item{x}{Numeric vector or a matrix}
	\item{eps}{Numerical. Shrinkage parameter of eigenvalue in 
	\code{ginverse_sym}}
	\item{a}{Numeric vector}
	\item{lambda}{Numeric value}
	\item{A}{Matrix}
    \item{y}{Matrix}
    \item{rho}{Numeric}
    \item{z}{Numeric}
}

%\details{
%xxx
%%}


\examples{
#############################################################################
## EXAMPLE 1: Trace of a matrix
#############################################################################	

set.seed(86)
A <- matrix( stats::runif(4) , 2 ,2 )
tracemat(A)
sum(diag(A))	# = sirt::tracemat(A)

#############################################################################
## EXAMPLE 2: Power function
#############################################################################

x <- 2.3
a <- 1.7
pow(x=x,a=a)
x^a			# = sirt::pow(x,a)

#############################################################################
## EXAMPLE 3: Soft and hard thresholding function (e.g. in LASSO estimation)
#############################################################################

x <- seq( -2 , 2 , length=100 )
y <- sirt::soft_thresholding( x , lambda = .5)
graphics::plot( x , y , type="l")

z <- sirt::hard_thresholding( x , lambda = .5)
graphics::lines( x , z , lty=2 , col=2)

#############################################################################
## EXAMPLE 4: Bounds on parameters
#############################################################################

pars <- c( .721 , .346 )
bounds_parameters( pars = pars , lower=c(-Inf, .5) , upper = c(Inf,1) )
}
