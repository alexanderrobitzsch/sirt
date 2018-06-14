%% File Name: stratified.cronbach.alpha.Rd
%% File Version: 0.19

\name{stratified.cronbach.alpha}
\alias{stratified.cronbach.alpha}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stratified Cronbach's Alpha
}
\description{
This function computes the stratified Cronbach's Alpha
for composite scales (Cronbach, Schoenemann & McKie, 1965;
Meyer, 2010).
}
\usage{
stratified.cronbach.alpha(data, itemstrata=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
An \eqn{N \times I} data frame
}
  \item{itemstrata}{
A matrix with two columns defining the item stratification.
The first column contains the item names, the second column
the item stratification label (these can be integers).
The default \code{NULL} does only compute Cronbach's Alpha
for the whole scale.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Cronbach, L.J., Schoenemann, P., & McKie, D. (1965).
Alpha coefficient for stratified-parallel tests.
\emph{Educational and Psychological Measurement},
\bold{25}, 291-312.

Meyer, P. (2010). \emph{Reliability}. Cambridge:
Oxford University Press.
}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
#############################################################################
# EXAMPLE 1: data.read
#############################################################################

data( data.read )
dat <- data.read
I <- ncol(dat)

# apply function without defining item strata
stratified.cronbach.alpha( data.read  )

# define item strata
itemstrata <- cbind( colnames(dat) , substring( colnames(dat) , 1 ,1 ) )
stratified.cronbach.alpha( data.read , itemstrata=itemstrata )
  ##   scale  I alpha mean.tot var.tot alpha.stratified
  ## 1 total 12 0.677    8.680   5.668            0.703
  ## 2     A  4 0.545    2.616   1.381               NA
  ## 3     B  4 0.381    2.811   1.059               NA
  ## 4     C  4 0.640    3.253   1.107               NA

\dontrun{
#**************************
# reliability analysis in psych package
library(psych)
# Cronbach's alpha and item discriminations
psych::alpha( dat )
# McDonald's omega
psych::omega(dat , nfactors=1)     # 1 factor
  ##   Alpha:                 0.69
  ##   Omega Total            0.69
##  => Note that alpha in this function is the standardized Cronbach's
##     alpha, i.e. alpha computed for standardized variables.
psych::omega(dat , nfactors=2)     # 2 factors
  ##   Omega Total            0.72
psych::omega(dat , nfactors=3)     # 3 factors
  ##   Omega Total            0.74
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Reliability}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
