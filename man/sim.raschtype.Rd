%% File Name: sim.raschtype.Rd
%% File Version: 0.16

\name{sim.raschtype}
\alias{sim.raschtype}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate from Generalized Logistic Item Response Model
}
\description{
This function simulates dichotomous item responses from a
generalized logistic item response model (Stukel, 1988).
The four-parameter logistic item response model
(Loken & Rulison, 2010) is a special case. See \code{\link{rasch.mml2}}
for more details.
}
\usage{
sim.raschtype(theta, b, alpha1 = 0, alpha2 = 0, fixed.a = NULL, 
    fixed.c = NULL, fixed.d = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
Unidimensional ability vector \eqn{\theta}
}
  \item{b}{
Vector of item difficulties \eqn{b}
}
  \item{alpha1}{
Parameter \eqn{\alpha_1} in generalized logistic link function
}
  \item{alpha2}{
Parameter \eqn{\alpha_2} in generalized logistic link function
}
  \item{fixed.a}{
Vector of item slopes \eqn{a}}
  \item{fixed.c}{
Vector of lower item asymptotes \eqn{c}
}
  \item{fixed.d}{
Vector of lower item asymptotes \eqn{d}
}
}
\details{
The class of generalized logistic link functions contain 
the most important link functions using the specifications (Stukel, 1988):

logistic link function: \eqn{\alpha_1=0} and \eqn{\alpha_2=0} \cr 
probit link function: \eqn{\alpha_1=0.165} and \eqn{\alpha_2=0.165} \cr
loglog link function: \eqn{\alpha_1=-0.037} and \eqn{\alpha_2=0.62} \cr
cloglog link function: \eqn{\alpha_1=0.62} and \eqn{\alpha_2=-0.037} \cr

See \code{\link{pgenlogis}} for exact transformation formulas of
the mentioned link functions.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Loken, E., & Rulison, K. L. (2010). Estimation of a four-parameter 
item response theory model. \emph{British Journal of Mathematical 
and Statistical Psychology}, \bold{63}, 509-525.

Stukel, T. A. (1988). Generalized logistic models. 
\emph{Journal of the American Statistical Association},
\bold{83}, 426-431.
}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{rasch.mml2}}, \code{\link{pgenlogis}}
}

\examples{
#############################################################################
## EXAMPLE 1: Simulation of data from a Rasch model (alpha_1 = alpha_2 = 0)
#############################################################################    

set.seed(9765)    
N <- 500    # number of persons
I <- 11     # number of items
b <- seq( -2 , 2 , length=I )
dat <- sirt::sim.raschtype( stats::rnorm( N ) , b )
colnames(dat) <- paste0( "I" , 1:I )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Simulating IRT models}
%\keyword{ ~kwd1 }

