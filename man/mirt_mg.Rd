%% File Name: mirt_mg.Rd
%% File Version: 0.037

\name{mirt_mg}
\alias{mirt_mg}
\alias{coef.mirt_mg}
\alias{logLik.mirt_mg}
\alias{vcov.mirt_mg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wrapper to \code{mirt} for Groupwise Fitting of the 2PL Model
}
\description{
Wrapper function to \code{mirt} for fitting the 2PL model separately for
groups.
}
\usage{
mirt_mg(resp, group, weights=NULL, ...)

\method{coef}{mirt_mg}(object,...)

\method{logLik}{mirt_mg}(object,...)

\method{vcov}{mirt_mg}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resp}{
Matrix of item responses
}
  \item{group}{
Group indicator variable
}
  \item{weights}{
Optional vector of sampling weights
}
  \item{\dots}{
Further arguments to be passed to \code{mirt}
}
\item{object}{Object of class \code{mirt_mg}}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
List with following entries
\item{coef}{Estimated item parameters}
\item{varmat}{Estimated variance matrix}
\item{infomat}{Estimated information matrix}
\item{parinfo}{Parameter information table}
\item{itempars}{Object which can be used in subsequent linking methods}
\item{logLik}{Sum of log-likelihood values of fitted objects}
\item{modlist}{List of fitted \code{mirt} objects}
\item{\dots}{Further elements}
}

%\references{
%% ~put references to the literature/web site here ~
%}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[mirt:mirt]{mirt::mirt}}
}

\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Fit 2PL model to the two gender groups in data.math
#############################################################################

data(data.math, package="sirt")

#- fit 2PL model in the two gender groups
rmod <- sirt::mirt_mg(resp=data.math$data[, -c(1,2)], group=data.math$data$female)
coef(rmod)
logLik(rmod)
vcov(rmod)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
