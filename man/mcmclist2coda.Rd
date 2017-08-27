%% File Name: mcmclist2coda.Rd
%% File Version: 0.10
%% File Last Change: 2017-01-18 18:08:39

\name{mcmclist2coda}
\alias{mcmclist2coda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write Coda File from an Object of Class \code{mcmc.list}
}
\description{
This function  writes a coda file from an object of class \code{mcmc.list}.
Note that only first entry (i.e. one chain) will be processed. 
}
\usage{
mcmclist2coda(mcmclist, name, coda.digits = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmclist}{
An object of class \code{mcmc.list}.
}
  \item{name}{
Name of the coda file to be written
}
  \item{coda.digits}{
Number of digits after decimal in the coda file
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The coda file and a corresponding index file are
written into the working directory.
}
%\references{
%% ~put references to the literature/web site here ~
%}
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
\dontrun{
#############################################################################
# EXAMPLE 1: MCMC estimation 2PNO dataset Reading
#############################################################################	

data(data.read)
# estimate 2PNO with MCMC with 3000 iterations and 500 burn-in iterations
mod <- sirt::mcmc.2pno( dat=data.read , iter=3000 , burnin=500 )
# plot MCMC chains
plot( mod$mcmcobj , ask=TRUE )
# write sampled chains into codafile
mcmclist2coda( mod$mcmcobj , name = "dataread_2pl" ) 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{coda}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
