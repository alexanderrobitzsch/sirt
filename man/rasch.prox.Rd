%% File Name: rasch.prox.Rd
%% File Version: 1.14

\name{rasch.prox}
\alias{rasch.prox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PROX Estimation Method for the Rasch Model
}
\description{
This function estimates the Rasch model using the PROX algorithm
(cited in Wright & Stone, 1999).
}
\usage{
rasch.prox(dat, dat.resp=1 - is.na(dat), freq=rep(1,nrow(dat)),
    conv=0.001, maxiter=30, progress=FALSE)
}
\arguments{
  \item{dat}{
An \eqn{N \times I} data frame of dichotomous response data. \code{NA}s
are not allowed and must be indicated by zero entries in the
response indicator matrix \code{dat.resp}.
}
  \item{dat.resp}{
An \eqn{N \times I} indicator data frame of nonmissing item responses.
}
  \item{freq}{
A vector of frequencies (or weights) of all rows in data frame \code{dat}.
}
  \item{conv}{
Convergence criterion for item parameters
}
  \item{maxiter}{
Maximum number of iterations
}
  \item{progress}{
Display progress?
}
}
%\details{
%%For details on the algorithm see Wright and Stone (1999, p. 135ff.).
%}
\value{
A list with following entries
\item{b}{Estimated item difficulties}
\item{theta}{Estimated person abilities}
\item{iter}{Number of iterations}
\item{sigma.i}{Item standard deviations}
\item{sigma.n}{Person standard deviations}
}
\references{
Wright, B., & Stone, W. (1999). \emph{Measurement Essentials}.
Wilmington: Wide Range.
}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
#############################################################################
# EXAMPLE 1: PROX data.read
#############################################################################

data(data.read)
mod <- sirt::rasch.prox( data.read )
mod$b       # item difficulties
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{PROX algorithm}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
