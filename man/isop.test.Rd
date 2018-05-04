%% File Name: isop.test.Rd
%% File Version: 0.15

\name{isop.test}
\alias{isop.test}
\alias{summary.isop.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Testing the ISOP Model
}
\description{
This function performs tests of the W1 axiom of the ISOP model
(Scheiblechner, 2003). Standard errors of the corresponding \eqn{W1_i} statistics
are obtained by Jackknife.
}

\usage{
isop.test(data, jackunits = 20, weights = rep(1, nrow(data)))

\method{summary}{isop.test}(object,...)   
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data frame with item responses
}
  \item{jackunits}{
A number of Jackknife units (if an integer is provided as the argument
value) or a vector in the Jackknife units are already
defined.
}
  \item{weights}{
Optional vector of sampling weights
}
\item{object}{Object of class \code{isop.test}}
\item{\dots}{Further arguments to be passed}
}



%\details{
%%  ~~ If necessary, more details than the description above ~~
%}



\value{
A list with following entries

  \item{itemstat}{Data frame with test and item statistics for the W1 axiom.
  The \eqn{W1_i} statistic is denoted as \code{est} while \code{se}
  is the corresponding standard error of the statistic.
  The sample size per item is \code{N} and \code{M} denotes the item mean.}
  \item{Es}{Number of concordances per item}
  \item{Ed}{Number of disconcordances per item}
  
The \eqn{W1_i} statistics are printed by the \code{summary} method.  
}

\references{
Scheiblechner, H. (2003). Nonparametric IRT: Testing the bi-isotonicity of
isotonic probabilistic models (ISOP). \emph{Psychometrika}, \bold{68},
79-96.
}

\author{
Alexander Robitzsch
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Fit the ISOP model with \code{\link{isop.dich}} or \code{\link{isop.poly}}.

See also the \pkg{ISOP} package at \emph{Rforge}: \emph{http://www.rforge.net/ISOP/}. 
}


\examples{
#############################################################################
# EXAMPLE 1: ISOP model data.Students
#############################################################################

data(data.Students, package="CDM")
dat <- data.Students[ , paste0("act",1:5) ]
dat <- dat[1:300 , ]    # select first 300 students

# perform the ISOP test
mod <- sirt::isop.test(dat)
summary(mod)
  ## -> W1i statistics
  ##     parm   N     M   est    se      t
  ##   1 test 300    NA 0.430 0.036 11.869
  ##   2 act1 278 0.601 0.451 0.048  9.384
  ##   3 act2 275 0.473 0.473 0.035 13.571
  ##   4 act3 274 0.277 0.352 0.098  3.596
  ##   5 act4 291 1.320 0.381 0.054  7.103
  ##   6 act5 276 0.460 0.475 0.042 11.184
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ISOP model}
\keyword{summary}
\keyword{Nonparametric IRT}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
