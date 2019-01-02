%% File Name: np.dich.Rd
%% File Version: 0.21

\name{np.dich}
\alias{np.dich}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nonparametric Estimation of Item Response Functions
}
\description{
This function does nonparametric item response function
estimation (Ramsay, 1991).
}
\usage{
np.dich(dat, theta, thetagrid, progress=FALSE, bwscale=1.1,
       method="normal")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
An \eqn{N \times I} data frame of dichotomous
item responses
}
  \item{theta}{
Estimated theta values, for example weighted likelihood
estimates from \code{\link{wle.rasch}}
}
  \item{thetagrid}{
A vector of theta values where the nonparametric item
response functions shall be evaluated.
}
  \item{progress}{
Display progress?
}
  \item{bwscale}{
The bandwidth parameter \eqn{h} is calculated by
the formula \eqn{h=}\code{bwscale}\eqn{\cdot N^{-1/5}}
}
  \item{method}{
The default \code{normal} performs kernel regression
with untransformed item responses. The method \code{binomial}
uses nonparametric logistic regression implemented
in the \pkg{sm} library.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with following entries
\item{dat}{Original data frame}
\item{thetagrid}{Vector of theta values at which the item response
    functions are evaluated}
\item{theta}{Used theta values as person parameter estimates}
\item{estimate}{Estimated item response functions}
\item{\dots}{}
}
\references{
Ramsay, J. O. (1991). Kernel smoothing approaches to
nonparametric item characteristic curve estimation.
\emph{Psychometrika, 56}, 611-630.
}

%\author{
%Alexander Robitzsch
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
# EXAMPLE 1: Reading dataset
#############################################################################
data( data.read )
dat <- data.read

# estimate Rasch model
mod <- sirt::rasch.mml2( dat )
# WLE estimation
wle1 <- sirt::wle.rasch( dat=dat, b=mod$item$b )$theta
# nonparametric function estimation
np1 <- sirt::np.dich( dat=dat, theta=wle1, thetagrid=seq(-2.5, 2.5, len=100 ) )
print( str(np1))
# plot nonparametric item response curves
plot( np1, b=mod$item$b )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Nonparametric item response theory}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
