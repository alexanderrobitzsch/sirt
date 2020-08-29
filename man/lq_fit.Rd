%% File Name: lq_fit.Rd
%% File Version: 0.16

\name{lq_fit}
\alias{lq_fit}
\alias{lq_fit_estimate_power}
\alias{dexppow}
\alias{rexppow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit of a \eqn{L_q} Regression Model
}
\description{
Fits a regression model in the \eqn{L_q} norm (also labeled as the \eqn{L_p} norm).
In more detail,
the optimization function \eqn{ \sum_i | y_i - x_i \beta | ^p} is optimized.
The nondifferentiable function is approximated by a differentiable approximation,
i.e., we use \eqn{|x| \approx \sqrt{x^2 + \varepsilon } }. The power \eqn{p}
can also be estimated by using \code{est_pow=TRUE}, see
Giacalone, Panarello and Mattera (2018). The algorithm iterates between estimating
regression coefficients and the estimation of power values. The estimation of the
power based on a vector of residuals \code{e} can be conducted using the
function \code{lq_fit_estimate_power}.

Using the \eqn{L_q} norm in the regression is equivalent to assuming an expontial
power function for residuals (Giacalone et al., 2018). The density function and
a simulation function is provided by \code{dexppow} and \code{rexppow}, respectively.
See also the \pkg{normalp} package.
}
\usage{
lq_fit(y, X, w=NULL, pow=2, eps=0.001, beta_init=NULL, est_pow=FALSE, optimizer="optim",
    eps_vec=10^seq(0,-10, by=-.5), conv=1e-4, miter=20, lower_pow=.1, upper_pow=5)

lq_fit_estimate_power(e, pow_init=2, lower_pow=.1, upper_pow=10)

dexppow(x, mu=0, sigmap=1, pow=2, log=FALSE)

rexppow(n, mu=0, sigmap=1, pow=2, xbound=100, xdiff=.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Dependent variable
}
  \item{X}{
Design matrix
}
  \item{w}{
Optional vector of weights
}
  \item{pow}{
Power \eqn{p} in \eqn{L_q} norm
}
  \item{est_pow}{
Logical indicating whether power should be estimated
}
  \item{eps}{
Parameter governing the differentiable approximation
}
\item{e}{Vector of resiuals}
\item{pow_init}{Initial value of power}
  \item{beta_init}{
Initial vector
}
  \item{optimizer}{
Can be \code{"optim"} or \code{"nlminb"}.
}
\item{eps_vec}{Vector with decreasing \eqn{\varepsilon} values used in
optimization}
\item{conv}{Convergence criterion}
\item{miter}{Maximum number of iterations}
\item{lower_pow}{Lower bound for estimated power}
\item{upper_pow}{Upper bound for estimated power}
\item{x}{Vector}
\item{mu}{Location parameter}
\item{sigmap}{Scale parameter}
\item{log}{Logical indicating whether the logarithm should be provided}
\item{n}{Sample size}
\item{xbound}{Lower and upper bound for density approximation}
\item{xdiff}{Grid width for density approximation}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
List with following several entries
\item{coefficients}{Vector of coefficients}
\item{res_optim}{Results of optimization}
\item{\ldots}{More values}
}

\references{
Giacalone, M., Panarello, D., & Mattera, R. (2018).
Multicollinearity in regression: an efficiency comparison between $L_p$-norm and least
squares estimators. \emph{Quality & Quantity, 52}(4), 1831-1859.
\doi{10.1007/s11135-017-0571-y}
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
#############################################################################
# EXAMPLE 1: Small simulated example with fixed power
#############################################################################

set.seed(98)
N <- 300
x1 <- stats::rnorm(N)
x2 <- stats::rnorm(N)
par1 <- c(1,.5,-.7)
y <- par1[1]+par1[2]*x1+par1[3]*x2 + stats::rnorm(N)
X <- cbind(1,x1,x2)

#- lm function in stats
mod1 <- stats::lm.fit(y=y, x=X)

#- use lq_fit function
mod2 <- sirt::lq_fit( y=y, X=X, pow=2, eps=1e-4)
mod1$coefficients
mod2$coefficients

\dontrun{
#############################################################################
# EXAMPLE 2: Example with estimated power values
#############################################################################

#*** simulate regression model with residuals from the exponential power distribution
#*** using a power of .30
set.seed(918)
N <- 2000
X <- cbind( 1, c(rep(1,N), rep(0,N)) )
e <- sirt::rexppow(n=2*N, pow=.3, xdiff=.01, xbound=200)
y <- X \%*\% c(1,.5) + e

#*** estimate model
mod <- sirt::lq_fit( y=y, X=X, est_pow=TRUE, lower_pow=.1)
mod1 <- stats::lm( y ~ 0 + X )
mod$coefficients
mod$pow
mod1$coefficients
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
