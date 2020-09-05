%% File Name: ccov.np.Rd
%% File Version: 0.232

\name{ccov.np}
\alias{ccov.np}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nonparametric Estimation of Conditional Covariances of Item Pairs
}
\description{
This function estimates conditional covariances of itempairs
(Stout, Habing, Douglas & Kim, 1996; Zhang & Stout,
1999a). The function is used for the estimation of the DETECT index.
The \code{ccov.np} function has the (default) option to smooth item response
functions (argument \code{smooth}) in the computation of conditional covariances
(Douglas, Kim, Habing, & Gao, 1998).
}
\usage{
ccov.np(data, score, bwscale=1.1, thetagrid=seq(-3, 3, len=200),
    progress=TRUE, scale_score=TRUE, adjust_thetagrid=TRUE, smooth=TRUE,
    use_sum_score=FALSE, bias_corr=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
An \eqn{N \times I} data frame of dichotomous responses.
Missing responses are allowed.
}
  \item{score}{
An ability estimate, e.g. the WLE
}
  \item{bwscale}{
Bandwidth factor for calculation of conditional covariance. The bandwidth
used in the estimation is \code{bwscale} times \eqn{N^{-1/5}}.
}
  \item{thetagrid}{
A vector which contains theta values where conditional
covariances are evaluated.
}
  \item{progress}{
Display progress?
}
\item{scale_score}{Logical indicating whether \code{score}
should be z standardized in advance of the calculation of
conditional covariances
}
\item{adjust_thetagrid}{Logical indicating whether \code{thetagrid} should be
adjusted if observed values in \code{score} are outside of \code{thetagrid}.
}
\item{smooth}{Logical indicating whether smoothing should be
applied for conditional covariance estimation}
\item{use_sum_score}{Logical indicating whether sum score should be used.
With this option, the bias corrected conditional covariance of Zhang and
Stout (1999) is used.}
\item{bias_corr}{Logical indicating whether bias correction (Zhang & Stout, 1999)
should be utilized if \code{use_sum_score=TRUE}.}
}

\references{
Douglas, J., Kim, H. R., Habing, B., & Gao, F. (1998). Investigating local dependence
with conditional covariance functions.
\emph{Journal of Educational and Behavioral Statistics, 23}(2), 129-151.
\doi{10.3102/10769986023002129}

Stout, W., Habing, B., Douglas, J., & Kim, H. R. (1996).
Conditional covariance-based nonparametric multidimensionality assessment.
\emph{Applied Psychological Measurement, 20}(4), 331-354.
\doi{10.1177/014662169602000403}

Zhang, J., & Stout, W. (1999). Conditional covariance structure
of generalized compensatory multidimensional items.
\emph{Psychometrika, 64}(2), 129-152.
\doi{10.1007/BF02294532}
}

%\author{
%Alexander Robitzsch
%}

\note{
This function is used in \code{\link{conf.detect}} and \code{\link{expl.detect}}.
}


\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: data.read | different settings for computing conditional covariance
#############################################################################

data(data.read, package="sirt")
dat <- data.read

#* fit Rasch model
mod <- sirt::rasch.mml2(dat)
score <- sirt::wle.rasch(dat=dat, b=mod$item$b)$theta

#* ccov with smoothing
cmod1 <- sirt::ccov.np(data=dat, score=score, bwscale=1.1)
#* ccov without smoothing
cmod2 <- sirt::ccov.np(data=dat, score=score, smooth=FALSE)

#- compare results
100*cbind( cmod1$ccov.table[1:6, "ccov"], cmod2$ccov.table[1:6, "ccov"])
}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{DETECT}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
