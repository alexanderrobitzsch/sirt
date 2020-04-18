%% File Name: lq_fit.Rd
%% File Version: 0.08

\name{lq_fit}
\alias{lq_fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit of a \eqn{L_q} Regression Model
}
\description{
Fits a regression model in the \eqn{L_q} norm. In more detail,
the optimization function \eqn{ \sum_i | y_i - x_i \beta | ^p} is optimized.
The nondifferentiable function is approximated by a differentiable approximation,
i.e., we use \eqn{|x| \approx \sqrt{x^2 + \varepsilon } }.
}
\usage{
lq_fit(y, X, w=NULL, pow=2, eps=0.001, beta_init=NULL, optimizer="optim",
    eps_vec=10^seq(0,-10, by=-.5))
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
  \item{eps}{
Parameter governing the differentiable approximation
}
  \item{beta_init}{
Initial vector
}
  \item{optimizer}{
Can be \code{"optim"} or \code{"nlminb"}.
}
\item{eps_vec}{Vector with decreasing \eqn{\varepsilon} values used in
optimization}
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

%\references{
%% ~put references to the literature/web site here ~
%}
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
# EXAMPLE 1: Small simulated example
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
