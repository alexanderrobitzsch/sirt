%% File Name: detect.index.Rd
%% File Version: 0.12
%% File Last Change: 2017-01-18 18:08:37

\name{detect.index}
\alias{detect.index}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of the DETECT and polyDETECT Index
}
\description{
This function calculated the DETECT and polyDETECT index (Stout, Habing, Douglas 
& Kim, 1996; Zhang & Stout, 1999a; Zhang, 2007). At first, 
conditional covariances have to be estimated
using the \code{\link{ccov.np}} function.
}
\usage{
detect.index(ccovtable, itemcluster)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ccovtable}{
A value of \code{\link{ccov.np}}.
}
  \item{itemcluster}{
Item cluster for each item. The order of entries must correspond
to the columns in \code{data} (submitted to \code{\link{ccov.np}}).
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Stout, W., Habing, B., Douglas, J., & Kim, H. R. (1996). 
Conditional covariance-based nonparametric multidimensionality assessment. 
\emph{Applied Psychological Measurement}, \bold{20}, 331-354.

Zhang, J., & Stout, W. (1999a). Conditional covariance structure 
of generalized compensatory multidimensional items. 
\emph{Psychometrika}, \bold{64}, 129-152.

Zhang, J., & Stout, W. (1999b). The theoretical DETECT index of 
dimensionality and its application to approximate simple structure. 
\emph{Psychometrika}, \bold{64}, 213-249. 

Zhang, J. (2007). Conditional covariance theory and DETECT for
polytomous items. \emph{Psychometrika}, \bold{72}, 69-91.
}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For examples see \code{\link{conf.detect}}.
}
%\examples{
%% ##
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{DETECT}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
