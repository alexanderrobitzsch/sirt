%% File Name: matrixfunctions.sirt.Rd
%% File Version: 1.161

\name{matrixfunctions.sirt}

\alias{rowMaxs.sirt}
\alias{rowMins.sirt}
\alias{rowCumsums.sirt}
\alias{colCumsums.sirt}
\alias{rowIntervalIndex.sirt}
\alias{rowKSmallest.sirt}
\alias{rowKSmallest2.sirt}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Some Matrix Functions
}
\description{
Some matrix functions which are written in \pkg{Rcpp} for speed
reasons.
}
\usage{
rowMaxs.sirt(matr)      # rowwise maximum

rowMins.sirt(matr)      # rowwise minimum

rowCumsums.sirt(matr)   # rowwise cumulative sum

colCumsums.sirt(matr)   # columnwise cumulative sum

rowIntervalIndex.sirt(matr,rn) # first index in row nn when matr(nn,zz) > rn(nn)

rowKSmallest.sirt(matr, K, break.ties=TRUE) # k smallest elements in a row
rowKSmallest2.sirt(matr, K )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matr}{
A numeric matrix
}
\item{rn}{A vector, usually a random number in applications}
\item{K}{An integer indicating the number of smallest elements to be
extracted}
\item{break.ties}{A logical which indicates if ties are randomly
    broken. The default is \code{TRUE}.
}
}

\details{
The function \code{rowIntervalIndex.sirt} searches for all rows \code{n}
the first index \code{i} for which \code{matr(n,i) > rn(n)} holds.

The functions \code{rowKSmallest.sirt} and \code{rowKSmallest2.sirt}
extract the \eqn{K} smallest entries in a matrix row. For small
numbers of \eqn{K} the function \code{rowKSmallest2.sirt} is
the faster one.
}

\value{
The output of \code{rowMaxs.sirt} is a list with the elements
\code{maxval} (rowwise maximum values) and \code{maxind} (rowwise
maximum indices). The output of \code{rowMins.sirt} contains
corresponding minimum values with entries \code{minval} and
\code{minind}.

The output of \code{rowKSmallest.sirt} are two matrices:
\code{smallval} contains the \eqn{K} smallest values whereas
\code{smallind} contains the \eqn{K} smallest indices.
}


%\references{
%% ~put references to the literature/web site here ~
%}


\author{
Alexander Robitzsch

The \pkg{Rcpp} code for \code{rowCumsums.sirt} is copied from code of
Romain Francois
(\url{https://lists.r-forge.r-project.org/pipermail/rcpp-devel/2010-October/001198.html}).
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For other matrix functions see the \pkg{matrixStats} package.
}


\examples{
#############################################################################
# EXAMPLE 1: a small toy example (I)
#############################################################################
set.seed(789)
N1 <- 10 ; N2 <- 4
M1 <- round( matrix( runif(N1*N2), nrow=N1, ncol=N2), 1 )

rowMaxs.sirt(M1)      # rowwise maximum
rowMins.sirt(M1)      # rowwise minimum
rowCumsums.sirt(M1)   # rowwise cumulative sum

# row index for exceeding a certain threshold value
matr <- M1
matr <- matr / rowSums( matr )
matr <- sirt::rowCumsums.sirt( matr )
rn <- runif(N1)    # generate random numbers
rowIntervalIndex.sirt(matr,rn)

# select the two smallest values
rowKSmallest.sirt(matr=M1, K=2)
rowKSmallest2.sirt(matr=M1, K=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Matrix utilities}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
