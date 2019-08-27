%% File Name: resp_groupwise.Rd
%% File Version: 0.06

\name{resp_groupwise}
\alias{resp_groupwise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates Group-Wise Item Response Dataset
}
\description{
Creates group-wise item response dataset.
}
\usage{
resp_groupwise(resp, group, items_group)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resp}{
Dataset with item responses
}
  \item{group}{
Vector of group identifiers
}
  \item{items_group}{
List containing vectors of groups for each item which should be made group-specific
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Dataset
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Toy dataset
#############################################################################

library(CDM)
library(TAM)

data(data.ex11, package="TAM")
dat <- data.ex11
dat[ dat==9 ] <- 0
resp <- dat[,-1]

# group labels
booklets <- sort( unique(paste(dat$booklet)))

#- fit initial model
mod0 <- TAM::tam.mml( resp, group=dat$booklet)
summary(mod0)

# fit statistics
fmod <- IRT.RMSD(mod)
stat <- abs(fmod$MD[,-1])
stat[ is.na( fmod$RMSD[,2:4] ) ] <- NA
thresh <- .01
round(stat,3)
# define list define groups for group-specific items
items_group <- apply( stat, 1, FUN=function(ll){
                v1 <- booklets[ which( ll > thresh ) ]
                v1[ ! is.na(v1) ]  } )

#- create extended response dataset
dat2 <- sirt::resp_groupwise(resp=resp, group=paste(dat$booklet), items_group=items_group)
colSums( ! is.na(dat2) )

#- fit model for extended response dataset
mod2 <- TAM::tam.mml( dat2, group=dat$booklet)
summary(mod2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{R utilities}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
