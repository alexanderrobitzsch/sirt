%% File Name: plot.mcmc.sirt.Rd
%% File Version: 0.15

\name{plot.mcmc.sirt}
\alias{plot.mcmc.sirt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Function for Objects of Class \code{mcmc.sirt}
}
\description{
Plot function for objects of class \code{mcmc.sirt}. These objects are generated
by: \code{\link{mcmc.2pno}}, \code{\link{mcmc.2pnoh}},
\code{\link{mcmc.3pno.testlet}}, \code{\link{mcmc.2pno.ml}}
}
\usage{
\method{plot}{mcmc.sirt}( x, layout=1, conflevel=0.9, round.summ=3,
   lag.max=.1, col.smooth="red", lwd.smooth=2, col.ci="orange",
   cex.summ=1, ask=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class \code{mcmc.sirt}
}
  \item{layout}{
Layout type. \code{layout=1} is the standard coda plot output,
\code{layout=2} gives a slightly different display.
}
  \item{conflevel}{
Confidence level (only applies to \code{layout=2})
}
  \item{round.summ}{
Number of digits to be rounded in summary (only applies to \code{layout=2})
}
  \item{lag.max}{
Maximum lag for autocorrelation plot (only applies to \code{layout=2}).
The default of .1 means that it is set to 1/10 of the number of iterations.
}
  \item{col.smooth}{
Color of smooth trend in traceplot (only applies to \code{layout=2})
}
  \item{lwd.smooth}{
Line type of smooth trend in traceplot (only applies to \code{layout=2})
}
  \item{col.ci}{
Color for displaying confidence interval (only applies to \code{layout=2})
}
  \item{cex.summ}{
Cex size for descriptive summary (only applies to \code{layout=2})
}
  \item{ask}{
Ask for a new plot (only applies to \code{layout=2})
}
  \item{\dots}{
Further arguments to be passed
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
\code{\link{mcmc.2pno}}, \code{\link{mcmc.2pnoh}},
\code{\link{mcmc.3pno.testlet}}, \code{\link{mcmc.2pno.ml}}
}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{plot}
% \keyword{methods}% __ONLY ONE__ keyword per line
