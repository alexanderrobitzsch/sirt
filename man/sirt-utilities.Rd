%% File Name: sirt-utilities.Rd
%% File Version: 0.408

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
\alias{sirt_sum_norm}
\alias{sirt_dnorm_discrete}
\alias{sirt_abs_smooth}
\alias{sirt_permutations}
\alias{sirt_rcpp_discrete_inverse}
\alias{move_variables_df}
\alias{print_digits}


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

#* normalize vector to have sum of one
sirt_sum_norm(x, na.rm=TRUE)
#* discrete normal distribution
sirt_dnorm_discrete(x, mean=0, sd=1, ...)

# plyr::rbind.fill implementation in sirt
sirt_rbind_fill(x, y)

# Fisher-z transformation, see psych::fisherz
sirt_fisherz(rho)
# inverse Fisher-z transformation, see psych::fisherz2r
sirt_antifisherz(z)

# smooth approximation of the absolute value function
sirt_abs_smooth(x, deriv=0, eps=1e-4)

# permutations with replacement
sirt_permutations(r,v)
  #-> is equivalent to gtools::permutations(n=length(v), r=D, v=v, repeats.allowed=TRUE)

# attach all elements in a list in a specified environment
sirt_attach_list_elements(x, envir)

# switch between stats::optim and stats::nlminb
sirt_optimizer(optimizer, par, fn, grad=NULL, method="L-BFGS-B", hessian=TRUE,
                   control=list(), ...)

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

# print a data frame x with fixed numbers of digits after the decimal
print_digits(x, digits=NULL)

# discrete inverse function
sirt_rcpp_discrete_inverse(x0, y0, y)

# move variables in a data frame
move_variables_df(x, after_var, move_vars)}

\arguments{
\item{pars}{Numeric vector}
\item{lower}{Numeric vector}
\item{upper}{Numeric vector}
\item{x}{Numeric vector or a matrix or a list}
\item{eps}{Numerical. Shrinkage parameter of eigenvalue in \code{ginverse_sym}}
\item{a}{Numeric vector}
\item{lambda}{Numeric value}
\item{A}{Matrix}
\item{nrow}{Integer}
\item{na.rm}{Logical}
\item{mean}{Numeric}
\item{sd}{Numeric}
\item{y}{Matrix}
\item{rho}{Numeric}
\item{deriv}{Integer indicating the order of derivative}
\item{z}{Numeric}
\item{r}{Integer}
\item{v}{Vector}
\item{envir}{Environment}
\item{optimizer}{Can be one of the following optimizers: \code{optim}, \code{nlminb},
\code{bobyqa} (from the \pkg{minqa} packlage), \code{Rvmmin} (from the
\pkg{optimx} package) or \code{nloptr} (from the \pkg{nloptr} package using
the argument \code{opts$algorithm="NLOPT_LD_MMA"}).
}
\item{par}{Initial parameter}
\item{fn}{Function}
\item{grad}{Gradient function}
\item{method}{Optimization method}
\item{hessian}{Logical}
\item{control}{Control list for R optimizers}
\item{\ldots}{Further arguments to be passed}
\item{obji}{Data frame}
\item{from}{Integer}
\item{to}{Integer}
\item{digits}{Integer}
\item{rownames_null}{Logical}
\item{grep_string}{String}
\item{pack}{Package name}
\item{CALL}{Call statement}
\item{x0}{Vector}
\item{y0}{Vector}
\item{after_var}{String indicating variable name after which variable specified
variables in \code{move_vars} should be moved}
\item{move_vars}{Variables which should be moved after \code{after_var}}
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

x <- seq(-2, 2, length=100)
y <- sirt::soft_thresholding( x, lambda=.5)
graphics::plot( x, y, type="l")

z <- sirt::hard_thresholding( x, lambda=.5)
graphics::lines( x, z, lty=2, col=2)

#############################################################################
## EXAMPLE 4: Bounds on parameters
#############################################################################

pars <- c(.721, .346)
bounds_parameters( pars=pars, lower=c(-Inf, .5), upper=c(Inf,1) )

#############################################################################
## EXAMPLE 5: Smooth approximation of absolute value function
#############################################################################

x <- seq(-1,1,len=100)
graphics::plot(x, abs(x), lwd=2, col=1, lty=1, type="l", ylim=c(-1,1) )
# smooth approximation
tt <- 2
graphics::lines(x, sirt::sirt_abs_smooth(x), lty=tt, col=tt, lwd=2)
# first derivative
tt <- 3
graphics::lines(x, sirt::sirt_abs_smooth(x, deriv=1), lty=tt, col=tt, lwd=2)
# second derivative
tt <- 4
graphics::lines(x, sirt::sirt_abs_smooth(x, deriv=2), lty=tt, col=tt, lwd=2)

# analytic computation of first and second derivative
stats::deriv( ~ sqrt(x^2 + eps), namevec="x", hessian=TRUE )

\dontrun{
#############################################################################
## EXAMPLE 6: Permutations with replacement
#############################################################################

D <- 4
v <- 0:1
sirt::sirt_permutations(r=D, v=v)
gtools::permutations(n=length(v), r=D, v=v, repeats.allowed=TRUE)
}
}
