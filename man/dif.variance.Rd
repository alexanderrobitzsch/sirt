\name{dif.variance}
\alias{dif.variance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DIF Variance
}
\description{
This function calculates the variance of DIF effects, the so called
DIF variance (Longford, Holland & Thayer, 1993).
}
\usage{
dif.variance(dif, se.dif, items = paste("item", 1:length(dif), sep = "") )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dif}{
Vector of uniform DIF effects
}
  \item{se.dif}{
Standard error of uniform DIF effects
}
  \item{items}{
Optional vector of item names
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with following entries
\item{weighted.DIFSD}{Weighted DIF standard deviation}
\item{unweigted.DIFSD}{DIF standard deviation}
\item{mean.se.dif}{Mean of standard errors of DIF effects}
\item{eb.dif}{Empirical Bayes estimates of DIF effects}
}
\references{
Longford, N. T., Holland, P. W., & Thayer, D. T. (1993).
Stability of the MH D-DIF statistics across populations.
In P. W. Holland & H. Wainer (Eds.). \emph{Differential
Item Functioning} (pp. 171-196). Hillsdale, NJ: Erlbaum.
}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{dif.logistic.regression}} for examples.
}
%\examples{
%% ### example
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Differential item functioning (DIF)}
\keyword{DIF variance}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
