%% File Name: md.pattern.sirt.Rd
%% File Version: 0.17

\name{md.pattern.sirt}
\alias{md.pattern.sirt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Response Pattern in a Binary Matrix
}
\description{
Computes different statistics of the response pattern in a binary
matrix.
}
\usage{
md.pattern.sirt(dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
A binary data matrix
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
A list with following entries:
\item{dat}{Original dataset}
\item{dat.resp1}{Indices for responses of 1's}
\item{dat.resp0}{Indices for responses of 0's}
\item{resp_patt}{Vector of response patterns}
\item{unique_resp_patt}{Unique response patterns}
\item{unique_resp_patt_freq}{Frequencies of unique response patterns}
\item{unique_resp_patt_firstobs}{First observation in original dataset
    \code{dat} of a unique response pattern}
\item{freq1}{Frequencies of 1's}
\item{freq0}{Frequencies of 0's}
\item{dat.ordered}{Dataset according to response patterns}
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

\seealso{
See also the \code{md.pattern} function in the \pkg{mice} package.
}

\examples{
#############################################################################
# EXAMPLE 1: Response patterns
#############################################################################
set.seed(7654)
N <- 21         # number of rows
I <- 4          # number of columns
dat <- matrix( 1*( stats::runif(N*I) > .3 ), N, I )
res <- sirt::md.pattern.sirt(dat)
# plot of response patterns
res$dat.ordered
image( z=t(res$dat.ordered), y=1:N, x=1:I, xlab="Items", ylab="Persons")
# 0's are yellow and 1's are red

#############################################################################
# EXAMPLE 2: Item response patterns for dataset data.read
#############################################################################

data(data.read)
dat <- data.read  ; N <- nrow(dat) ; I <- ncol(dat)
# order items according to p values
dat <- dat[, order(colMeans(dat, na.rm=TRUE )) ]
# analyzing response pattern
res <- sirt::md.pattern.sirt(dat)
res$dat.ordered
image( z=t(res$dat.ordered), y=1:N, x=1:I, xlab="Items", ylab="Persons")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Utilities}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
