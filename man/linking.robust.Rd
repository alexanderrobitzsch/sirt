%% File Name: linking.robust.Rd
%% File Version: 0.25

\name{linking.robust}
\alias{linking.robust}
\alias{summary.linking.robust}
\alias{plot.linking.robust}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robust Linking of Item Intercepts
}

\description{
This function implements a robust alternative of mean-mean linking which
employs trimmed means instead of means.
The linking constant is calculated for varying trimming parameters \eqn{k}.
}

\usage{
linking.robust(itempars)

\method{summary}{linking.robust}(object,...)

\method{plot}{linking.robust}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{itempars}{
Data frame of item parameters (item intercepts). The first column
contains the item label, the 2nd and 3rd columns item parameters of
two studies.
}
\item{object}{Object of class \code{linking.robust}}
\item{x}{Object of class \code{linking.robust}}
\item{\dots}{Further arguments to be passed}
}



\value{
A list with following entries

\item{ind.kopt}{Index for optimal scale parameter}
\item{kopt}{Optimal scale parameter}
\item{meanpars.kopt}{Linking constant for optimal scale parameter}
\item{se.kopt}{Standard error for linking constant obtained with optimal scale parameter}
\item{meanpars}{Linking constant dependent on the scale parameter}
\item{se}{Standard error of the linking constant dependent on the scale parameter}
\item{sd}{DIF standard deviation (non-robust estimate)}
\item{mad}{DIF standard deviation (robust estimate using the MAD measure)}
\item{pars}{Original item parameters}
\item{k.robust}{Used vector of scale parameters}
\item{I}{Number of items}
\item{itempars}{Used data frame of item parameters}
}

%\references{
%% ~put references to the literature/web site here ~
%}


%\author{
%Alexander Robitzsch
%}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\seealso{
Other functions for linking: \code{\link{linking.haberman}},
\code{\link{equating.rasch}}

See also the \pkg{plink} package.
}
%% ~~objects to See Also as \code{\link{help}}, ~~~

\examples{
#############################################################################
# EXAMPLE 1: Linking data.si03
#############################################################################

data(data.si03)
res1 <- sirt::linking.robust( itempars=data.si03 )
summary(res1)
  ##   Number of items=27
  ##   Optimal trimming parameter k=8 |  non-robust parameter k=0
  ##   Linking constant=-0.0345 |  non-robust estimate=-0.056
  ##   Standard error=0.0186 |  non-robust estimate=0.027
  ##   DIF SD: MAD=0.0771 (robust) | SD=0.1405 (non-robust)
plot(res1)

\dontrun{
#############################################################################
# EXAMPLE 2: Linking PISA item parameters data.pisaPars
#############################################################################

data(data.pisaPars)

# Linking with items
res2 <- sirt::linking.robust( data.pisaPars[, c(1,3,4)] )
summary(res2)
  ##   Optimal trimming parameter k=0 |  non-robust parameter k=0
  ##   Linking constant=-0.0883 |  non-robust estimate=-0.0883
  ##   Standard error=0.0297 |  non-robust estimate=0.0297
  ##   DIF SD: MAD=0.1824 (robust) | SD=0.1487 (non-robust)
##  -> no trimming is necessary for reducing the standard error
plot(res2)

#############################################################################
# EXAMPLE 3: Linking with simulated item parameters containing outliers
#############################################################################

# simulate some parameters
I <- 38
set.seed(18785)
itempars <- data.frame("item"=paste0("I",1:I) )
itempars$study1 <- stats::rnorm( I, mean=.3, sd=1.4 )
# simulate DIF effects plus some outliers
bdif <- stats::rnorm(I,mean=.4,sd=.09)+( stats::runif(I)>.9 )* rep( 1*c(-1,1)+.4, each=I/2 )
itempars$study2 <- itempars$study1 + bdif

# robust linking
res <- sirt::linking.robust( itempars )
summary(res)
  ##   Number of items=38
  ##   Optimal trimming parameter k=12 |  non-robust parameter k=0
  ##   Linking constant=-0.4285 |  non-robust estimate=-0.5727
  ##   Standard error=0.0218 |  non-robust estimate=0.0913
  ##   DIF SD: MAD=0.1186 (robust) | SD=0.5628 (non-robust)
## -> substantial differences of estimated linking constants in this case of
##    deviations from normality of item parameters
plot(res)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Linking}
%% \keyword{Robust linking}% __ONLY ONE__ keyword per line
%% \keyword{summary}
%% \keyword{plot}
