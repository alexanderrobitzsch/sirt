%% File Name: equating.rasch.jackknife.Rd
%% File Version: 0.24

\name{equating.rasch.jackknife}
\alias{equating.rasch.jackknife}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Jackknife Equating Error in Generalized Logistic Rasch Model
}
\description{
This function estimates the linking error in linking
based on Jackknife (Monseur & Berezner, 2007).
}
\usage{
equating.rasch.jackknife(pars.data, display=TRUE,
   se.linkerror=FALSE, alpha1=0, alpha2=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pars.data}{
Data frame with four columns: jackknife unit (1st column),
item parameter study 1 (2nd column), item parameter study 2 (3rd column),
item (4th column)
}
  \item{display}{
Display progress?
}
  \item{se.linkerror}{
Compute standard error of the linking error
}
  \item{alpha1}{
Fixed \eqn{\alpha_1} parameter in the generalized item
response model
}
  \item{alpha2}{
Fixed \eqn{\alpha_2} parameter in the generalized item
response model
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with following entries:
\item{pars.data}{Used item parameters}
\item{itemunits}{Used units for jackknife}
\item{descriptives}{Descriptives for Jackknife.
\code{linkingerror.jackknife} is the estimated linking error.}
}
\references{
Monseur, C., & Berezner, A. (2007). The computation of equating errors
in international surveys in education.
\emph{Journal of Applied Measurement, 8}, 323-335.
}

%\author{
%Alexander Robitzsch
%}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For more details on linking methods see \code{\link{equating.rasch}}.
}
\examples{
#############################################################################
# EXAMPLE 1: Linking errors PISA study
#############################################################################

data(data.pisaPars)
pars <- data.pisaPars

# Linking error: Jackknife unit is the testlet
vars <- c("testlet","study1","study2","item")
res1 <- sirt::equating.rasch.jackknife(pars[, vars])
res1$descriptives
  ##   N.items N.units      shift        SD linkerror.jackknife SE.SD.jackknife
  ## 1      25       8 0.09292838 0.1487387          0.04491197      0.03466309

# Linking error: Jackknife unit is the item
res2 <- sirt::equating.rasch.jackknife(pars[, vars ] )
res2$descriptives
  ##   N.items N.units      shift        SD linkerror.jackknife SE.SD.jackknife
  ## 1      25      25 0.09292838 0.1487387          0.02682839      0.02533327
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Equating}
%% \keyword{Linking}% __ONLY ONE__ keyword per line
