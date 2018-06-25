%% File Name: gom.jml.Rd
%% File Version: 0.20

\name{gom.jml}
\alias{gom.jml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Grade of Membership Model (Joint Maximum Likelihood Estimation)
}
\description{
This function estimates the grade of membership model
employing a joint maximum likelihood estimation
method (Erosheva, 2002; p. 23ff.).
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\usage{
gom.jml(dat, K=2, seed=NULL, globconv=0.001, maxdevchange=0.001,
        maxiter=600, min.lambda=0.001, min.g=0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
    Data frame of dichotomous item responses
}
  \item{K}{
    Number of classes
}
\item{seed}{
    Seed value of random number generator. Deterministic starting values
    are used for the default value \code{NULL}.
}
  \item{globconv}{
    Global parameter convergence criterion
}
  \item{maxdevchange}{
    Maximum change in relative deviance
}
  \item{maxiter}{
    Maximum number of iterations
}
  \item{min.lambda}{
    Minimum \eqn{\lambda_{ik}} parameter to be estimated
}
  \item{min.g}{
    Minimum \eqn{g_{pk}} parameter to be estimated
}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETAILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\details{
The item response model of the grade of membership model with \eqn{K} classes
for dichotomous correct responses \eqn{X_{pi}}
of person \eqn{p} on item \eqn{i} is
\deqn{
    P(X_{pi}=1 | g_{p1}, \ldots, g_{pK} )=\sum_k \lambda_{ik} g_{pk}
\quad, \quad \sum_k g_{pk}=1
                }{%
    P(X_{pi}=1 | g_{p1}, \ldots, g_{pK} )=\sum_k \lambda_{ik} g_{pk}
,    \sum_k g_{pk}=1
            }
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\value{
A list with following entries:
\item{lambda}{Data frame of item parameters \eqn{\lambda_{ik}}}
\item{g}{Data frame of individual membership scores \eqn{g_{pk}}}
\item{g.mean}{Mean membership scores}
\item{gcut}{Discretized membership scores}
\item{gcut.distr}{Distribution of discretized membership scores}
\item{K}{Number of classes}
\item{deviance}{Deviance}
\item{ic}{Information criteria}
\item{N}{Number of students}
\item{score}{Person score}
\item{iter}{Number of iterations}
\item{datproc}{List with processed data (recoded data, starting values, ...)}
\item{\dots}{Further values}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



\references{
Erosheva, E. A. (2002). \emph{Grade of membership
and latent structure models with application to
disability survey data}.
PhD thesis, Carnegie Mellon University, Department
of Statistics.
}

\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
S3 method \code{\link{summary.gom}}
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\examples{
#############################################################################
# EXAMPLE 1: TIMSS data
#############################################################################

data( data.timss)
dat <- data.timss$data[, grep("M", colnames(data.timss$data) ) ]

# 2 Classes (deterministic starting values)
m2 <- sirt::gom.jml(dat,K=2, maxiter=10 )
summary(m2)

\dontrun{
# 3 Classes with fixed seed and maximum number of iterations
m3 <- sirt::gom.jml(dat,K=3, maxiter=50,seed=89)
summary(m3)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Grade of membership model}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
