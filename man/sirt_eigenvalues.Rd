%% File Name: sirt_eigenvalues.Rd
%% File Version: 0.15

\name{sirt_eigenvalues}
\alias{sirt_eigenvalues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
First Eigenvalues of a Symmetric Matrix
}
\description{
This function computes the first \eqn{D} eigenvalues and eigenvectors of a
symmetric positive definite matrices. The eigenvalues are computed
by the Rayleigh quotient method (Lange, 2010, p. 120).
}
\usage{
sirt_eigenvalues( X, D, maxit=200, conv=10^(-6) )
}

%% eigenvalues.sirt <- function (X,D,maxit,conv){
%%    .Call("eigenvalues_sirt", X,D,maxit,conv, PACKAGE="sirt")
%%                    }

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Symmetric matrix
}
  \item{D}{
Number of eigenvalues to be estimated
  }
  \item{maxit}{
Maximum number of iterations
}
  \item{conv}{
Convergence criterion
}
}

%% \details{
%% xxx
%% }

\value{
A list with following entries:
\item{d}{Vector of eigenvalues
}
\item{u}{Matrix with eigenvectors in columns}
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




\examples{
Sigma <- diag(1,3)
Sigma[ lower.tri(Sigma) ] <- Sigma[ upper.tri(Sigma) ] <- c(.4,.6,.8 )
sirt::sirt_eigenvalues(X=Sigma, D=2 )
# compare with svd function
svd(Sigma)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Eigenvalues}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
