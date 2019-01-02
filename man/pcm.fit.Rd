%% File Name: pcm.fit.Rd
%% File Version: 0.20

\name{pcm.fit}
\alias{pcm.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Item and Person Fit Statistics for the Partial Credit Model
}

\description{
Computes item and person fit statistics in the
partial credit model (Wright & Masters, 1990).
The rating scale model is accommodated as a particular partial
credit model (see Example 3).
}

\usage{
pcm.fit(b, theta, dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{
Matrix with item category parameters (see Examples)
}
  \item{theta}{
Vector with estimated person parameters
}
  \item{dat}{
Dataset with item responses
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
A list with entries
 \item{itemfit}{Item fit statistics}
 \item{personfit}{Person fit statistics}
}

\references{
Wright, B. D., & Masters, G. N. (1990). Computation of outfit and
infit statistics. \emph{Rasch Measurement Transactions, 3}:4, 84-85.
}

%\author{
%%  ~~who you are~~
%}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{personfit.stat} for person fit statistics for dichotomous
item responses. See also the \pkg{PerFit} package for further person
fit statistics.

Item fit in other \R packages:
\code{eRm::itemfit},
\code{\link[TAM:tam.fit]{TAM::tam.fit}},
\code{\link[mirt:itemfit]{mirt::itemfit}},
\code{ltm::item.fit},

Person fit in other \R packages:
\code{eRm::itemfit},
\code{\link[mirt:itemfit]{mirt::itemfit}},
\code{ltm::person.fit},

See \code{\link{pcm.conversion}} for conversions of different
parametrizations of the partial credit model.
}



\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Partial credit model
#############################################################################

data(data.Students,package="CDM")
dat <- data.Students
# select items
items <- c(paste0("sc", 1:4 ), paste0("mj", 1:4 ) )
dat <- dat[,items]
dat <- dat[ rowSums( 1 - is.na(dat) ) > 0, ]

#*** Model 1a: Partial credit model in TAM
# estimate model
mod1a <- TAM::tam.mml( resp=dat )
summary(mod1a)
# estimate person parameters
wle1a <- TAM::tam.wle(mod1a)
# extract item parameters
b1 <- - mod1a$AXsi[, -1 ]
# parametrization in xsi parameters
b2 <- matrix( mod1a$xsi$xsi, ncol=3, byrow=TRUE )
# convert b2 to b1
b1b <- 0*b1
b1b[,1] <- b2[,1]
b1b[,2] <- rowSums( b2[,1:2] )
b1b[,3] <- rowSums( b2[,1:3] )
# assess fit
fit1a <- sirt::pcm.fit(b=b1, theta=wle1a$theta, dat)
fit1a$item

#############################################################################
# EXAMPLE 2: Rasch model
#############################################################################

data(data.read)
dat <- data.read

#*** Rasch model in TAM
# estimate model
mod <- TAM::tam.mml( resp=dat )
summary(mod)
# estimate person parameters
wle <- TAM::tam.wle(mod)
# extract item parameters
b1 <- - mod$AXsi[, -1 ]
# assess fit
fit1a <- sirt::pcm.fit(b=b1, theta=wle$theta, dat)
fit1a$item

#############################################################################
# EXAMPLE 3: Rating scale model
#############################################################################

data(data.Students,package="CDM")
dat <- data.Students
items <- paste0("sc", 1:4 )
dat <- dat[,items]
dat <- dat[ rowSums( 1 - is.na(dat) ) > 0, ]

#*** Model 1: Rating scale model in TAM
# estimate model
mod1 <- tam.mml( resp=dat, irtmodel="RSM")
summary(mod1)
# estimate person parameters
wle1 <- tam.wle(mod1)
# extract item parameters
b1 <- - mod1a$AXsi[, -1 ]
# fit statistic
pcm.fit(b=b1, theta=wle1$theta, dat)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Item fit}
%% \keyword{Person fit}% __ONLY ONE__ keyword per line
