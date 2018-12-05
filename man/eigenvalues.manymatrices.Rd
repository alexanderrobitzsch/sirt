%% File Name: eigenvalues.manymatrices.Rd
%% File Version: 0.16

\name{eigenvalues.manymatrices}
\alias{eigenvalues.manymatrices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computation of Eigenvalues of Many Symmetric Matrices
}
\description{
This function computes the eigenvalue decomposition of \eqn{N}
symmetric positive definite matrices. The eigenvalues are computed
by the Rayleigh quotient method (Lange, 2010, p. 120). In addition,
the inverse matrix can be calculated.
}
\usage{
eigenvalues.manymatrices(Sigma.all, itermax=10, maxconv=0.001,
    inverse=FALSE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Sigma.all}{
An \eqn{N \times D^2} matrix containing the \eqn{D^2} entries
of \eqn{N} symmetric matrices of dimension \eqn{D \times D}
}
  \item{itermax}{
Maximum number of iterations
}
  \item{maxconv}{
Convergence criterion for convergence of eigenvectors
}
\item{inverse}{A logical which indicates if the inverse matrix
    shall be calculated}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with following entries
  \item{lambda}{Matrix with eigenvalues}
  \item{U}{An \eqn{N \times D^2} Matrix of orthonormal eigenvectors}
  \item{logdet}{Vector of logarithm of determinants}
  \item{det}{Vector of determinants}
  \item{Sigma.inv}{Inverse matrix if \code{inverse=TRUE}.}
%% ...
}
\references{
Lange, K. (2010). \emph{Numerical Analysis for Statisticians}.
New York: Springer.
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
# define matrices
Sigma <- diag(1,3)
Sigma[ lower.tri(Sigma) ] <- Sigma[ upper.tri(Sigma) ] <- c(.4,.6,.8 )
Sigma1 <- Sigma

Sigma <- diag(1,3)
Sigma[ lower.tri(Sigma) ] <- Sigma[ upper.tri(Sigma) ] <- c(.2,.1,.99 )
Sigma2 <- Sigma

# collect matrices in a "super-matrix"
Sigma.all <- rbind( matrix( Sigma1, nrow=1, byrow=TRUE),
                matrix( Sigma2, nrow=1, byrow=TRUE) )
Sigma.all <- Sigma.all[ c(1,1,2,2,1 ), ]

# eigenvalue decomposition
m1 <- sirt::eigenvalues.manymatrices( Sigma.all )
m1

# eigenvalue decomposition for Sigma1
s1 <- svd(Sigma1)
s1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Eigenvalues}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
