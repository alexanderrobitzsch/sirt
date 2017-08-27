%% File Name: monoreg.rowwise.Rd
%% File Version: 0.10
%% File Last Change: 2017-01-18 18:08:39

\name{monoreg.rowwise}
\alias{monoreg.rowwise}
\alias{monoreg.colwise}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Monotone Regression for Rows or Columns in a Matrix
}

\description{
Monotone (isotone) regression for rows (\code{monoreg.rowwise}) or
columns (\code{monoreg.colwise}) in a matrix.
}

\usage{
monoreg.rowwise(yM, wM)

monoreg.colwise(yM, wM)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{yM}{
Matrix with dependent variable for the regression.
Values are assumed to be sorted.
}
\item{wM}{
Matrix with weights for every entry in the \code{yM}
matrix.
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
Matrix with fitted values
}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Alexander Robitzsch

The \code{monoreg} function from the \pkg{fdrtool}
package is simply extended to handle matrix input.
}

\note{
This function is used for fitting the ISOP model 
(see \code{\link{isop.dich}}).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also the \code{monoreg} function from the \pkg{fdrtool}
package.
}


\examples{
y <- c(22.5 , 23.33 , 20.83 , 24.25 )
w <- c( 3,3,3,2)
# define matrix input
yM <- matrix( 0 , nrow=2 , ncol=4 )
wM <- yM
yM[1,] <- yM[2,] <- y
wM[1,] <- w
wM[2,] <- c(1,3,4, 3 )

# fit rowwise monotone regression
monoreg.rowwise( yM , wM )
# compare results with monoreg function from fdrtool package
\dontrun{
miceadds::library_install("fdrtool")
fdrtool::monoreg(x=yM[1,] , w=wM[1,])$yf
fdrtool::monoreg(x=yM[2,] , w=wM[2,])$yf
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Monotone regression}
\keyword{Isotone regression}
% __ONLY ONE__ keyword per line
