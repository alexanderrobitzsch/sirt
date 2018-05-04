%% File Name: R2noharm.jackknife.Rd
%% File Version: 1.11

\name{R2noharm.jackknife}
\alias{R2noharm.jackknife}
\alias{summary.R2noharm.jackknife}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Jackknife Estimation of NOHARM Analysis
}
\description{
This function performs a jackknife estimation of NOHARM analysis
to get standard errors based on a replication method (see Christoffersson, 1977).
}
\usage{
R2noharm.jackknife(object, jackunits = NULL)

\method{summary}{R2noharm.jackknife}(object, logfile=NULL , \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Object of class \code{R2noharm}
}
  \item{jackunits}{
A vector of integers or a number. If it is a number, then it refers
to the number of jackknife units. If it is a vector of integers, then this vector
defines the allocation of persons jackknife units. Integers corresponds to
row indexes in the data set.
}
\item{logfile}{File name if the summary should be sunk into a file}
\item{\dots}{Further arguments to be passed}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
A list of lists with following entries:

\item{partable}{
    Data frame with parameters
        }
\item{se.pars}{List of estimated standard errors for all parameter estimates:
    \code{tanaka.stat}, \code{rmsr.stat}, \code{rmsea.stat},
    \code{chisquare_df.stat}, \code{thresholds.stat}, \code{final.constants.stat},
    \code{uniquenesses.stat}, \code{factor.cor.stat}, \code{loadings.stat},
    \code{loadings.theta.stat}
    }
\item{jackknife.pars}{List with obtained results by jackknifing for all parameters:
    \code{j.tanaka}, \code{j.rmsr}, \code{rmsea}, \code{chisquare_df}, 
    \code{j.pm}, \code{j.thresholds}, \code{j.factor.cor}, 
    \code{j.loadings}, \code{j.loadings.theta}
        }
\item{u.jacknunits}{Unique jackknife elements}
}


\references{
Christoffersson, A. (1977). Two-step weighted least squares factor analysis of 
dichotomized variables. \emph{Psychometrika}, \bold{42}, 433-438.
}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{R2noharm}}
}
%\examples{
%# See examples for function R2noharm
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{NOHARM}
\keyword{summary}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
