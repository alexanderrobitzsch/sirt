%% File Name: mgsem.Rd
%% File Version: 0.082

\name{mgsem}
\alias{mgsem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of Multiple-Group Structural Equation Models
}
\description{
Estimates a multiple-group structural equation model. The function allows arbitrary
prior distributions on model parameters and regularized estimation with the SCAD and
the LASSO penalty. Moreover, it can also conduct robust moment estimation using
the \eqn{L_p} loss function \eqn{\rho(x)=|x|^p} for \eqn{p \ge 0}.
See Robitzsch (2023) for more details.
}
\usage{
mgsem(suffstat, model, data=NULL, group=NULL, weights=NULL, estimator="ML",
     p_me=2, p_pen=1, pen_type="scad", diffpar_pen=NULL, pen_sample_size=TRUE,
     a_scad=3.7, eps_approx=0.001, comp_se=TRUE, se_delta_formula=FALSE,
     prior_list=NULL, hessian=TRUE, fixed_parms=FALSE, cd=FALSE,
     cd_control=list(maxiter=20, tol=5*1e-04, interval_length=0.05, method="exact"),
     partable_start=NULL, num_approx=FALSE, technical=NULL, control=list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{suffstat}{
List containing sufficient statistics
}
  \item{model}{
Model specification, see examples. Can have entries \code{est}, \code{index},
\code{lower}, \code{upper}, \code{prior}, \code{pen_l2}, \code{pen_lp},
\code{pen_difflp}. Each entry can be defined for model matrices \code{ALPHA},
\code{NU}, \code{LAM}, \code{PHI}, and \code{PSI}.
}
  \item{data}{
Optional data frame
}
  \item{group}{
Optional vector of group identifiers
}
  \item{weights}{
Optional vector of sampling weights
}
  \item{estimator}{
Character. Can be either \code{"ML"} for maximum likelihood fitting function or
\code{"ME"} for robust moment estimation.
}
  \item{p_me}{
Power in $L_p$ loss function for robust moment estimation
}
  \item{p_pen}{
Power for penalty in regularized estimation. For regular LASSO and SCAD penalty
functions, it is $p=1$.
}
  \item{pen_type}{
Penalty type. Can be either \code{"scad"} or \code{"lasso"}.
}
  \item{diffpar_pen}{
List containing values of regularization parameters in fused lasso estimation
}
  \item{pen_sample_size}{
List containing values for sample sizes for regularized estimation
}
  \item{a_scad}{
Parameter $a$ used in SCAD penalty
}
  \item{eps_approx}{
Approximation value for nondifferentiable robust moment fitting function or
penalty function
}
  \item{comp_se}{
Logical indicating whether standard errors should be computed
}
  \item{se_delta_formula}{
Logical indicating whether standard errors should be computed according to the
delta formula
}
  \item{prior_list}{
List containing specifications of the prior distributions
}
  \item{hessian}{
Logical indicating whether the Hessian matrix should be computed
}
  \item{fixed_parms}{
Logical indicating whether all model parameters should be fixed
}
  \item{cd}{
Logical indicating whether coordinate descent should be used for estimation
}
  \item{cd_control}{
Control parameters for coordinate descent estimation
}
  \item{partable_start}{
Starting values for parameter estimation
}
  \item{num_approx}{
Logical indicating whether derivatives should be computed based on numerical
differentiation
}
  \item{technical}{
Parameters used for optimization in \code{sirt_optimizer}
}
  \item{control}{
Control paramaters for optimization
}
}

\details{
[MORE INFORMATION TO BE ADDED]
}


\value{
A list with following values
\item{coef}{Coeffients}
\item{vcov}{Variance matrix}
\item{se}{Vector of standard errors}
\item{partable}{Parameter table}
\item{model}{Specified model}
\item{opt_res}{Result from optimization}
\item{opt_value}{Value of fitting function}
\item{residuals}{Residuals of sufficient statistics}
\item{ic}{Information criteria}
\item{technical}{Specifications of optimizer}
\item{suffstat_vcov}{Variance matrix of sufficient statistics}
\item{data_proc}{Processed data}
\item{case_ll}{Casewise log-likelihood function}
\item{\dots}{Further values}
}

\references{
Robitzsch, A. (2023). Model-robust estimation of multiple-group structural
equation models. \emph{Algorithms, 16}(4), 210. \doi{10.3390/a16040210}
}


%\author{
%%  ~~who you are~~
%}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}


\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Noninvariant item intercepts in a multiple-group SEM
#############################################################################

#---- simulate data
set.seed(65)
G <- 3  # number of groups
I <- 5  # number of items
# define lambda and nu parameters
lambda <- matrix(1, nrow=G, ncol=I)
nu <- matrix(0, nrow=G, ncol=I)
err_var <- matrix(1, nrow=G, ncol=I)

# define extent of noninvariance
dif_int <- .5

#- 1st group: N(0,1)
nu[1,4] <- dif_int
#- 2nd group: N(0.3,1.5)
gg <- 2 ;
nu[gg,1] <- -dif_int
#- 3nd group: N(.8,1.2)
gg <- 3
nu[gg,2] <- -dif_int
#- define distributions of groups
mu <- c(0,.3,.8)
sigma <- sqrt(c(1,1.5,1.2))
N <- rep(1000,3) # sample sizes per group

exact <- FALSE
suffstat <- sirt::invariance_alignment_simulate(nu, lambda, err_var, mu, sigma, N,
                output="suffstat", groupwise=TRUE, exact=exact)

#---- model specification

# model specifications joint group
est <- list(
        ALPHA=matrix( c(0), ncol=1),
        NU=matrix( 0, nrow=I, ncol=1),
        LAM=matrix(1, nrow=I, ncol=1),
        PHI=matrix(0,nrow=1,ncol=1),
        PSI=diag(rep(1,I))
        )

# parameter index
index <- list(
        ALPHA=0*est$ALPHA,
        NU=1+0*est$NU,
        LAM=1+0*est$LAM,
        PHI=0*est$PHI,
        PSI=diag(1,I)
        )

# lower bounds
lower <- list(
        PSI=diag(rep(0.01,I)), PHI=matrix(0.01,1,1)
        )

#*** joint parameters
group0 <- list(est=est, index=index, lower=lower)

#*** group1
est <- list(
        ALPHA=matrix( c(0), ncol=1),
        NU=matrix( 0, nrow=I, ncol=1),
        LAM=matrix(0, nrow=I, ncol=1),
        PHI=matrix(1,nrow=1,ncol=1)
            )

# parameter index
index <- list(
        ALPHA=0*est$ALPHA,
        NU=0*est$NU,
        LAM=1*est$LAM,
        PHI=0*est$PHI
        )

group1 <- list(est=est, index=index, lower=lower)

#*** group 2 and group 3

# modify parameter index
index$ALPHA <- 1+0*est$ALPHA
index$PHI <- 1+0*est$PHI
group3 <- group2 <- list(est=est, index=index, lower=lower)

#*** define model
model <- list(group0=group0, group1=group1, group2=group2, group3=group3)

#-- estimate model with ML
res1 <- sirt::mgsem( suffstat=suffstat, model=model2, eps_approx=1e-4, estimator="ML",
                    technical=list(maxiter=500, optimizer="optim"),
                    hessian=FALSE, comp_se=FALSE, control=list(trace=1) )
str(res1)

#-- robust moment estimation with p=0.5

optimizer <- "optim"
technical <- list(maxiter=500, optimizer=optimizer)
eps_approx <- 1e-3

res2 <- sirt::mgsem( suffstat=suffstat, model=res1$model, p_me=0.5,
                    eps_approx=eps_approx, estimator="ME", technical=technical,
                    hessian=FALSE, comp_se=FALSE, control=list(trace=1) )

#---- regularized estimation

nu_lam <- 0.1    # regularization parameter

# redefine model
define_model <- function(model, nu_lam)
{
    pen_lp <- list( NU=nu_lam+0*model$group1$est$NU)
    ee <- "group1"
    for (ee in c("group1","group2","group3"))
    {
        model[[ee]]$index$NU <- 1+0*index$NU
        model[[ee]]$pen_lp <- pen_lp
    }
    return(model)
}

model3 <- define_model(model=model, nu_lam=nu_lam)
pen_type <- "scad"

res3 <- sirt::mgsem( suffstat=suffstat, model=model3, p_pen=1, pen_type=pen_type,
                    eps_approx=eps_approx, estimator="ML",
                    technical=list(maxiter=500, optimizer="optim"),
                    hessian=FALSE, comp_se=FALSE, control=list(trace=1) )
str(res3)
}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
