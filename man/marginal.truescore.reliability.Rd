\name{marginal.truescore.reliability}
\alias{marginal.truescore.reliability}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
True-Score Reliability for Dichotomous Data
}
\description{
This function computes the marginal true-score reliability for
dichotomous data (Dimitrov, 2003; May & Nicewander, 1994) for 
the four-parameter logistic item response model
(see \code{\link{rasch.mml2}} for details regarding this IRT model).
}
\usage{
marginal.truescore.reliability(b, a=1+0*b ,c=0*b ,d=1+0*b,
    mean.trait=0, sd.trait=1, theta.k=seq(-6,6,len=200) )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{
Vector of item difficulties
}
  \item{a}{
Vector of item discriminations
}
  \item{c}{
Vector of guessing parameters
}
  \item{d}{
Vector of upper asymptotes
}
  \item{mean.trait}{
Mean of trait distribution
}
  \item{sd.trait}{
Standard deviation of trait distribution
}
  \item{theta.k}{
Grid at which the trait distribution should be evaluated
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with following entries:
  \item{rel.test}{Reliability of the test}
  \item{item}{True score variance (\code{sig2.true}, error variance
  (\code{sig2.error}) and item reliability (\code{rel.item}).
  Expected proportions correct are in the column \code{pi}.
    }
  \item{pi}{Average proportion correct for all items and persons}
  \item{sig2.tau}{True score variance \eqn{\sigma^2_{\tau}}
  (calculated by the formula in May & Nicewander, 1994)}
  \item{sig2.error}{Error variance \eqn{\sigma^2_{e}}}
}
\references{
Dimitrov, D. (2003). Marginal true-score measures and reliability
for binary items as a function of their IRT parameters.
\emph{Applied Psychological Measurement}, \bold{27}, 440-458.

May, K., & Nicewander, W. A. (1994). Reliability and information
functions for percentile ranks. \emph{Journal of Educational
Mesaurement}, \bold{31}, 313-325.
}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{greenyang.reliability}} for calculating the reliability
for multidimensional measures.
}
\examples{
#############################################################################
# EXAMPLE 1: Dimitrov (2003) Table 1 - 2PL model
#############################################################################

# item discriminations
a <- 1.7*c(0.449,0.402,0.232,0.240,0.610,0.551,0.371,0.321,0.403,0.434,0.459,
    0.410,0.302,0.343,0.225,0.215,0.487,0.608,0.341,0.465)
# item difficulties
b <- c( -2.554,-2.161,-1.551,-1.226,-0.127,-0.855,-0.568,-0.277,-0.017,
    0.294,0.532,0.773,1.004,1.250,1.562,1.385,2.312,2.650,2.712,3.000 )

marginal.truescore.reliability( b=b , a =a )
  ##   Reliability= 0.606 

#############################################################################
# EXAMPLE 2: Dimitrov (2003) Table 2
#  3PL model: Poetry items (4 items)
#############################################################################

# slopes, difficulties and guessing parameters
a <- 1.7*c(1.169,0.724,0.554,0.706 )
b <- c(0.468,-1.541,-0.042,0.698 )
c <- c(0.159,0.211,0.197,0.177 )

res <- sirt::marginal.truescore.reliability( b=b , a =a , c=c)
  ##   Reliability= 0.403 
  ##   > round( res$item , 3 )
  ##     item    pi sig2.tau sig2.error rel.item
  ##   1    1 0.463    0.063      0.186    0.252
  ##   2    2 0.855    0.017      0.107    0.135
  ##   3    3 0.605    0.026      0.213    0.107
  ##   4    4 0.459    0.032      0.216    0.130

#############################################################################
# EXAMPLE 3: Reading Data
#############################################################################
data( data.read)

#***
# Model 1: 1PL
mod <- sirt::rasch.mml2( data.read )
marginal.truescore.reliability( b=mod$item$b )
  ##   Reliability= 0.653 

#***
# Model 2: 2PL
mod <- sirt::rasch.mml2( data.read , est.a=1:12 )
marginal.truescore.reliability( b=mod$item$b , a=mod$item$a)
  ##   Reliability= 0.696 

\dontrun{
# compare results with Cronbach's alpha and McDonald's omega
# posing a 'wrong model' for normally distributed data
library(psych)
psych::omega(dat , nfactors=1)     # 1 factor
  ##  Omega_h for 1 factor is not meaningful, just omega_t
  ##   Omega 
  ##   Call: omega(m = dat, nfactors = 1)
  ##   Alpha:                 0.69 
  ##   G.6:                   0.7 
  ##   Omega Hierarchical:    0.66 
  ##   Omega H asymptotic:    0.95 
  ##   Omega Total            0.69 

##! Note that alpha in psych is the standardized one. 
} 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Reliability}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
