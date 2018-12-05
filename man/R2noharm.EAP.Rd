%% File Name: R2noharm.EAP.Rd
%% File Version: 0.15

\name{R2noharm.EAP}
\alias{R2noharm.EAP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
EAP Factor Score Estimation
}
\description{
This function performs EAP factor score estimation of
an item response model estimated with NOHARM.
}
\usage{
R2noharm.EAP(noharmobj, theta.k=seq(-6, 6, len=21), print.output=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{noharmobj}{
Object of class \code{R2noharm} or \code{noharm.sirt}
}
  \item{theta.k}{
Vector of discretized theta values on which the posterior
is evaluated. This vector applies to all dimensions.
}
\item{print.output}{An optional logical indicating whether output should
be displayed at the console}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with following entries
  \item{person}{Data frame of person parameter EAP estimates and
  their corresponding standard errors
    }
  \item{theta}{Grid of multidimensional theta values where the posterior is
      evaluated.}
  \item{posterior}{Individual posterior distribution evaluated at \code{theta}}
  \item{like}{Individual likelihood}
  \item{EAP.rel}{EAP reliabilities of all dimensions}
  \item{probs}{Item response probabilities evaluated at \code{theta}}
}
%\references{
%% ~put references to the literature/web site here ~
%}


%\author{
%Alexander Robitzsch
%}



%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For examples see \code{\link{R2noharm}} and \code{\link{noharm.sirt}}.
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Factor scores}
\keyword{NOHARM}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
