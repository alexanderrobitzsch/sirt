%% File Name: isop.scoring.Rd
%% File Version: 1.26

\name{isop.scoring}
\alias{isop.scoring}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scoring Persons and Items in the ISOP Model
}
\description{
This function does the scoring in the isotonic probabilistic
model (Scheiblechner, 1995, 2003, 2007).
Person parameters are ordinally scaled but the ISOP
model also allows \emph{specific objective} (ordinal) comparisons
for persons (Scheiblechner, 1995).
}
\usage{
isop.scoring(dat,score.itemcat=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame with dichotomous or polytomous item responses
}
\item{score.itemcat}{
Optional data frame with scoring points for every item and
every category (see Example 2).
}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETAILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\details{
This function extracts the scoring rule of the ISOP model
(if \code{score.itemcat !=NULL}) and calculates the
modified percentile score for every person. The score
\eqn{s_{ik}} for item \eqn{i} and category \eqn{k}
is calculated as
\deqn{ s_{ik}=\sum_{j=0}^{k-1} f_{ij} -
\sum_{j=k+1}^K f_{ij}=P( X_i < k  )
- P( X_i > k ) }
where \eqn{f_{ik}} is the relative frequency of item \eqn{i}
in category \eqn{k} and \eqn{K} is the maximum category.
The modified percentile score \eqn{\rho_p}
for subject \eqn{p} (\code{mpsc} in \code{person}) is
defined by
\deqn{ \rho_p=\frac{1}{I} \sum_{i=1}^I
        \sum_{j=0}^K s_{ik} \mathbf{1}( X_{pi}=k ) }
Note that for dichotomous items, the sum score is a
sufficient statistic for \eqn{\rho_p} but this is
not the case for polytomous items.
The modified percentile score \eqn{\rho_p}
ranges between -1 and 1.

The modified item P-score \eqn{\rho_i} (Scheiblechner, 2007, p. 52) is
defined by
\deqn{ \rho_i=\frac{1}{I-1} \cdot \sum_j \left[ P( X_j < X_i )
        - P( X_j > X_i ) \right ] }
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\value{
A list with following entries:
  \item{person}{
  A data frame with person parameters. The modified
  percentile score \eqn{\rho_p} is denoted by \code{mpsc}.
      }
  \item{item}{Item statistics and scoring parameters.
    The item P-scores \eqn{\rho_i}
   are labeled as \code{pscore}.}
  \item{p.itemcat}{Frequencies for every item category}
  \item{score.itemcat}{
        Scoring points for every item category
                }
  \item{distr.fct}{Empirical distribution function}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\references{
Scheiblechner, H. (1995). Isotonic ordinal
probabilistic models (ISOP). \emph{Psychometrika,
60}, 281-304.

Scheiblechner, H. (2003). \emph{Nonparametric IRT:
Scoring functions and ordinal parameter estimation
of isotonic probabilistic models (ISOP)}.
Technical Report, Philipps-Universitaet Marburg.

Scheiblechner, H. (2007). A unified nonparametric IRT model for
d-dimensional psychological test data (d-ISOP).
\emph{Psychometrika, 72}, 43-67.
}



%\author{
%Alexander Robitzsch
%}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For fitting the ISOP and ADISOP model see
\code{\link{isop.dich}} or \code{\link{fit.isop}}.
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\examples{
#############################################################################
# EXAMPLE 1: Dataset Reading
#############################################################################

data( data.read )
dat <- data.read

# Scoring according to the ISOP model
msc <- sirt::isop.scoring( dat )
# plot student scores
boxplot( msc$person$mpsc ~ msc$person$score )

#############################################################################
# EXAMPLE 2: Dataset students from CDM package | polytomous items
#############################################################################

library("CDM")
data( data.Students, package="CDM")
dat <- stats::na.omit(data.Students[, -c(1:2) ])

# Scoring according to the ISOP model
msc <- sirt::isop.scoring( dat )
# plot student scores
boxplot( msc$person$mpsc ~ msc$person$score )

# scoring with known scoring rule for activity items
items <- paste0( "act", 1:5 )
score.itemcat <- msc$score.itemcat
score.itemcat <- score.itemcat[ items, ]
msc2 <- sirt::isop.scoring( dat[,items], score.itemcat=score.itemcat )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ISOP model}
%% \keyword{Nonparametric IRT}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
