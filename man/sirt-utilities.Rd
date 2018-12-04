%% File Name: sirt-utilities.Rd
%% File Version: 0.34

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
\alias{sirt_attach_list_elements}
\alias{sirt_summary_print_objects}
\alias{sirt_summary_print_package_rsession}
\alias{sirt_summary_print_package}
\alias{sirt_summary_print_rsession}
\alias{sirt_summary_print_call}
\alias{sirt_optimizer}
\alias{sirt_matrix2}
\alias{sirt_colMeans}
\alias{sirt_colSDs}
\alias{sirt_colMins}
\alias{sirt_colMaxs}
\alias{sirt_colMedians}

\title{Utility Functions in \pkg{sirt}}

\description{
Utility functions in \pkg{sirt}.
}

\usage{
# bounds entries in a vector
bounds_parameters( pars, lower=NULL, upper=NULL)

# improper density function which always returns a value of 1
dimproper(x)

# generalized inverse of a symmetric function
ginverse_sym(A, eps=1E-8)
# hard thresholding function
hard_thresholding(x, lambda)
# soft thresholding function
soft_thresholding(x, lambda)

# power function x^a, like in Cpp
pow(x, a)
# trace of a matrix
tracemat(A)

#** matrix functions
sirt_matrix2(x, nrow)   # matrix() function with byrow=TRUE
sirt_colMeans(x, na.rm=TRUE)
sirt_colSDs(x, na.rm=TRUE)
sirt_colMins(x, na.rm=TRUE)
sirt_colMaxs(x, na.rm=TRUE)
sirt_colMedians(x, na.rm=TRUE)

# plyr::rbind.fill implementation in sirt
sirt_rbind_fill(x, y)

# Fisher-z transformation, see psych::fisherz
sirt_fisherz(rho)
# inverse Fisher-z transformation, see psych::fisherz2r
sirt_antifisherz(z)

# attach all elements in a list in a specified environment
sirt_attach_list_elements(x, envir)

# switch between stats::optim and stats::nlminb
sirt_optimizer(optimizer, par, fn, grad=NULL, method="L-BFGS-B", hessian=TRUE, ...)

# print objects in a summary
sirt_summary_print_objects(obji, from=NULL, to=NULL, digits=3, rownames_null=TRUE,
      grep_string=NULL)
# print package version and R session
sirt_summary_print_package_rsession(pack)
# print package version
sirt_summary_print_package(pack)
# print R session
sirt_summary_print_rsession()
# print call
sirt_summary_print_call(CALL)
}

\arguments{
\item{pars}{Numeric vector}
\item{lower}{Numeric vector}
\item{upper}{Numeric vector}
\item{x}{Numeric vector or a matrix or a list}
\item{eps}{Numerical. Shrinkage parameter of eigenvalue in
\code{ginverse_sym}}
\item{a}{Numeric vector}
\item{lambda}{Numeric value}
\item{A}{Matrix}
\item{nrow}{Integer}    
\item{na.rm}{Logical}
\item{y}{Matrix}
\item{rho}{Numeric}
\item{z}{Numeric}
\item{envir}{Environment}
\item{optimizer}{Can be \code{optim} or \code{nlminb}}
\item{par}{Initial parameter}
\item{fn}{Function}
\item{grad}{Gradient function}
\item{method}{Optimization method}
\item{hessian}{Logical}
\item{\ldots}{Further arguments to be passed}
\item{obji}{Data frame}
\item{from}{Integer}
\item{to}{Integer}
\item{digits}{Integer}
\item{rownames_null}{Logical}
\item{grep_string}{String}
\item{pack}{Package name}
\item{CALL}{Call statement}
}

%\details{
%xxx
%%}


\examples{
#############################################################################
## EXAMPLE 1: Trace of a matrix
#############################################################################

set.seed(86)
A <- matrix( stats::runif(4), 2,2 )
tracemat(A)
sum(diag(A))    #=sirt::tracemat(A)

#############################################################################
## EXAMPLE 2: Power function
#############################################################################

x <- 2.3
a <- 1.7
pow(x=x,a=a)
x^a            #=sirt::pow(x,a)

#############################################################################
## EXAMPLE 3: Soft and hard thresholding function (e.g. in LASSO estimation)
#############################################################################

x <- seq( -2, 2, length=100 )
y <- sirt::soft_thresholding( x, lambda=.5)
graphics::plot( x, y, type="l")

z <- sirt::hard_thresholding( x, lambda=.5)
graphics::lines( x, z, lty=2, col=2)

#############################################################################
## EXAMPLE 4: Bounds on parameters
#############################################################################

pars <- c( .721, .346 )
bounds_parameters( pars=pars, lower=c(-Inf, .5), upper=c(Inf,1) )
}
