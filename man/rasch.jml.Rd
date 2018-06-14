%% File Name: rasch.jml.Rd
%% File Version: 1.32

\name{rasch.jml}
\alias{rasch.jml}
\alias{summary.rasch.jml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Joint Maximum Likelihood (JML) Estimation of the Rasch Model
}
\description{
This function estimates the Rasch model using joint maximum likelihood
estimation (Lincare, 1994). The PROX algorithm (Lincare, 1994) is used
for the generation of starting values of item parameters.
}
\usage{
rasch.jml(dat, method = "MLE", b.init = NULL, constraints = NULL, weights = NULL,
    glob.conv = 10^(-6),   conv1 = 1e-05, conv2 = 0.001, progress = TRUE,
    bsteps = 4,thetasteps = 2, wle.adj = 0, jmliter = 100, prox = TRUE ,
    proxiter = 30, proxconv = 0.01, dp=NULL , theta.init = NULL , calc.fit=TRUE)

\method{summary}{rasch.jml}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{An \eqn{N \times I} data frame of dichotomous item responses
      where \eqn{N} indicates the number of persons
and \eqn{I} the number of items
}
  \item{method}{
Method for estimating person parameters during JML iterations.
\code{MLE} is maximum likelihood estimation (where person with
perfect scores are deleted from analysis). \code{WLE} uses weighted
likelihood estimation (Warm, 1989) for person parameter estimation.
Default is \code{MLE}.
}
  \item{b.init}{
Initial values of item difficulties
}
  \item{constraints}{
Optional matrix or data.frame with two columns. First column is an integer of
item indexes or item names (\code{colnames(dat)}) which shall be fixed
during estimation. The second column is the corresponding
item difficulty.
}
  \item{weights}{
Person sample weights. Default is \code{NULL}, i.e. all persons in the sample
are equally weighted.
}
  \item{glob.conv}{
Global convergence criterion with respect to the log-likelihood function
}
  \item{conv1}{
Convergence criterion for estimation of item parameters
}
  \item{conv2}{
Convergence criterion for estimation of person parameters
}
  \item{progress}{
Display progress? Default is \code{TRUE}
}
%% bsteps = 4,thetasteps = 2, wle.adj = 0,
  \item{bsteps}{
Number of steps for b parameter estimation
}
  \item{thetasteps}{
Number of steps for theta parameter estimation
}
  \item{wle.adj}{
Score adjustment for WLE estimation
}
  \item{jmliter}{
Number of maximal iterations during JML estimation
}
  \item{prox}{
Should the PROX algorithm (see \code{\link{rasch.prox}}) be used as initial estimations?
Default is \code{TRUE}.
}
  \item{proxiter}{
Number of maximal PROX iterations
}
  \item{proxconv}{
Convergence criterion for PROX iterations
}
  \item{dp}{
 Object created from data preparation function (\code{.data.prep})
 which could be created in earlier JML runs. Default is \code{NULL}.
 }
\item{theta.init}{
 Initial person parameter estimate
 }
\item{calc.fit}{
 Should itemfit being calculated?
 }
  \item{object}{
Object of class \code{rasch.jml}
}
  \item{\dots}{
Further arguments to be passed
}
}

\details{
The estimation is known to have a bias in item parameters for
a fixed (finite) number of items. In literature (Lincare, 1994), a simple
bias correction formula is proposed and included in the value
\code{item$itemdiff.correction} in this function. If \eqn{I} denotes the number
of items, then the correction factor is \eqn{\frac{I-1}{I}}.
}
%%If there are \eqn{J}
%%item parameter constraints, then the correction factor amounts to
%%\eqn{\frac{I-J-1}{I-J}}.



\value{
A list with following entries
\item{item}{Estimated item parameters}
\item{person}{Estimated person parameters}
\item{method}{Person parameter estimation method}
\item{dat}{Original data frame}
\item{deviance}{Deviance}
\item{data.proc}{Processed data frames excluding persons with extreme scores}
\item{dp}{Value of data preparation (it is used in the function
\code{\link{rasch.jml.jackknife1}})}
}
\references{
 Linacre, J. M. (1994). \emph{Many-Facet Rasch Measurement}. Chicago: MESA Press.

 Warm, T. A. (1989). Weighted likelihood estimation of ability in the item
 response theory. \emph{Psychometrika}, \bold{54}, 427-450.
    }
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Get a summary with \code{\link{summary.rasch.jml}}.

See \code{\link{rasch.prox}} for the PROX algorithm as initial iterations.

For a bias correction of the JML method try \code{\link{rasch.jml.jackknife1}}.

JML estimation can also be conducted with the \pkg{TAM}
(\code{\link[TAM:tam.jml2]{TAM::tam.jml2}})
and \pkg{mixRasch} (\code{mixRasch::mixRasch})
packages.

See also marginal maximum likelihood estimation with \code{\link{rasch.mml2}}
or the \R package \pkg{ltm}.
}

% (\code{\link[mixRasch:mixRasch]{mixRasch::mixRasch}})


\examples{
#############################################################################
# EXAMPLE 1: Simulated data from the Rasch model
#############################################################################
set.seed(789)
N <- 500    # number of persons
I <- 11     # number of items
b <- seq( -2 , 2 , length=I )
dat <- sirt::sim.raschtype( stats::rnorm( N ) , b )
colnames(dat) <- paste( "I" , 1:I , sep="")

# JML estimation of the Rasch model
mod1 <- sirt::rasch.jml( dat )
summary(mod1)

# MML estimation with rasch.mml2 function
mod2 <- sirt::rasch.mml2( dat )
summary(mod2)

# Pairwise method of Fischer
mod3 <- sirt::rasch.pairwise( dat )
summary(mod3)

# JML estimation in TAM
\dontrun{
library(TAM)
mod4 <- TAM::tam.jml2( resp=dat )

# JML estimation in mixRasch package
library(mixRasch)
mod5 <- mixRasch::mixRasch( dat, steps=1, n.c=1, max.iter=50)
print(mod5)
mod5$item.par       # extract item parameters

#******
# item parameter constraints in JML estimation
# fix item difficulties: b[4]=-.76 and b[6]= .10
constraints <- matrix( cbind( 4 , -.76 ,
                              6 , .10 ) ,
                  ncol=2 , byrow=TRUE )
mod6 <- sirt::rasch.jml( dat , constraints = constraints )
summary(mod6)
  # For constrained item parameters, it this not obvious
  # how to calculate a 'right correction' of item parameter bias
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{Rasch model}
\keyword{Joint maximum likelihood (JML)}
\keyword{summary}
% __ONLY ONE__ keyword per line
