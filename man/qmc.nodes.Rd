%% File Name: qmc.nodes.Rd
%% File Version: 0.11

\name{qmc.nodes}
\alias{qmc.nodes}
\title{
Calculation of Quasi Monte Carlo Integration Points
}
\description{
This function calculates integration nodes based on the multivariate
normal distribution with zero mean vector and identity covariance
matrix. See Pan and Thompson (2007) and Gonzales et al. (2006)
for details.
}
\usage{
qmc.nodes(snodes, ndim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snodes}{
Number of integration nodes
}
  \item{ndim}{
Number of dimensions
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
 \item{theta}{A matrix of integration points}
}
\references{
Gonzalez, J., Tuerlinckx, F., De Boeck, P., & Cools, R. (2006). 
Numerical integration in logistic-normal models. 
\emph{Computational Statistics & Data Analysis}, \bold{51}, 1535-1548.	
	
Pan, J., & Thompson, R. (2007). Quasi-Monte Carlo estimation in
generalized linear mixed models. \emph{Computational Statistics &
Data Analysis}, \bold{51}, 5765-5775.
}
\author{
Alexander Robitzsch
}

\note{
This function uses the
\code{\link[sfsmisc:QUnif]{sfsmisc::QUnif}} function from
the \pkg{sfsmisc} package.
}
% sfsmisc::QUnif

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
## some toy examples

# 5 nodes on one dimension
qmc.nodes( snodes=5 , ndim=1 )
  ##            [,1]
  ## [1,]  0.0000000
  ## [2,] -0.3863753
  ## [3,]  0.8409238
  ## [4,] -0.8426682
  ## [5,]  0.3850568

# 7 nodes on two dimensions
qmc.nodes( snodes =7 , ndim=2 )
  ##             [,1]        [,2]
  ## [1,]  0.00000000 -0.43072730
  ## [2,] -0.38637529  0.79736332
  ## [3,]  0.84092380 -1.73230641
  ## [4,] -0.84266815 -0.03840544
  ## [5,]  0.38505683  1.51466109
  ## [6,] -0.00122394 -0.86704605
  ## [7,]  1.35539115  0.33491073
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Quasi Monte Carlo integration}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
