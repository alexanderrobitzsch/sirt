%% File Name: personfit.stat.Rd
%% File Version: 0.23

\name{personfit.stat}
\alias{personfit.stat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Person Fit Statistics for the Rasch Model
}
\description{
This function collects some person fit statistics
for the Rasch model (Karabatsos, 2003; Meijer & Sijtsma, 2001).
}
\usage{
personfit.stat(dat, abil, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
An \eqn{N \times I} data frame of dichotomous item
responses
}
  \item{abil}{
An ability estimate, e.g. the WLE
}
  \item{b}{
Estimated item difficulty
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A data frame with following columns (see Meijer & Sijtsma 2001
for a review of different person fit statistics):
\item{case}{Case index}
\item{abil}{Ability estimate \code{abil}}
\item{mean}{Person mean of correctly solved items}
\item{caution}{Caution index}
\item{depend}{Dependability index}
\item{ECI1}{\eqn{ECI1}}
\item{ECI2}{\eqn{ECI2}}
\item{ECI3}{\eqn{ECI3}}
\item{ECI4}{\eqn{ECI4}}
\item{ECI5}{\eqn{ECI5}}
\item{ECI6}{\eqn{ECI6}}
\item{l0}{Fit statistic \eqn{l_0}}
\item{lz}{Fit statistic \eqn{l_z}}
\item{outfit}{Person outfit statistic}
\item{infit}{Person infit statistic}
\item{rpbis}{Point biserial correlation of item responses
and item \eqn{p} values}
\item{rpbis.itemdiff}{Point biserial correlation of item responses
and item difficulties \code{b}}
\item{U3}{Fit statistic \eqn{U_3}}
}

\references{
Karabatsos, G. (2003). Comparing the aberrant response detection performance
of thirty-six person-fit statistics. \emph{Applied Measurement in Education,
16}, 277-298.

Meijer, R. R., & Sijtsma, K. (2001). Methodology
review: Evaluating person fit. \emph{Applied Psychological
Measurement, 25}, 107-135.
}

\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{pcm.fit}} for person fit in the partial credit model.

See the \pkg{irtProb} and \pkg{PerFit} packages for person fit statistics
and person response curves and functions included in other packages:
\code{\link[mirt:personfit]{mirt::personfit}},
\code{eRm::personfit} and
\code{ltm::person.fit}.
}

% \code{\link[TAM:tam.fa]{tam.fa}} (\pkg{TAM})

\examples{
#############################################################################
# EXAMPLE 1: Person fit Reading Data
#############################################################################

data(data.read)
dat <- data.read

# estimate Rasch model
mod <- sirt::rasch.mml2( dat )
# WLE
wle1 <- sirt::wle.rasch( dat,b=mod$item$b )$theta
b <- mod$item$b # item difficulty

# evaluate person fit
pf1 <- sirt::personfit.stat( dat=dat, abil=wle1, b=b)

\dontrun{
# dimensional analysis of person fit statistics
x0 <- stats::na.omit(pf1[, -c(1:3) ] )
stats::factanal( x=x0, factors=2, rotation="promax" )
  ## Loadings:
  ##                Factor1 Factor2
  ## caution         0.914
  ## depend          0.293   0.750
  ## ECI1            0.869   0.160
  ## ECI2            0.869   0.162
  ## ECI3            1.011
  ## ECI4            1.159  -0.269
  ## ECI5            1.012
  ## ECI6            0.879   0.130
  ## l0              0.409  -1.255
  ## lz             -0.504  -0.529
  ## outfit          0.297   0.702
  ## infit           0.362   0.695
  ## rpbis          -1.014
  ## rpbis.itemdiff  1.032
  ## U3              0.735   0.309
  ##
  ## Factor Correlations:
  ##         Factor1 Factor2
  ## Factor1   1.000  -0.727
  ## Factor2  -0.727   1.000
  ##
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Person fit}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
