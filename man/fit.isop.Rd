%% File Name: fit.isop.Rd
%% File Version: 1.20

\name{fit.isop}
\alias{fit.isop}
\alias{fit.adisop}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Fitting the ISOP and ADISOP Model for Frequency Tables
}
\description{
Fit the isotonic probabilistic model (ISOP;
Scheiblechner, 1995) and the additive isotonic probabilistic
model (ADISOP; Scheiblechner, 1999).
}

\usage{
fit.isop(freq.correct, wgt, conv=1e-04, maxit=100,
      progress=TRUE, calc.ll=TRUE)

fit.adisop(freq.correct, wgt, conv=1e-04, maxit=100,
      epsilon=0.01, progress=TRUE, calc.ll=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{freq.correct}{
Frequency table
}
  \item{wgt}{
Weights for frequency table (number of persons in
each cell)
}
  \item{conv}{
Convergence criterion
}
  \item{maxit}{
Maximum number of iterations
}
\item{epsilon}{Additive constant to handle cell frequencies
of 0 or 1 in \code{fit.adisop}
}
  \item{progress}{
Display progress?
}
\item{calc.ll}{Calculate log-likelihood values?
The default is \code{TRUE}.}
}
\details{
See \code{\link{isop.dich}} for more details of the
ISOP and ADISOP model.
}


\value{
A list with following entries
\item{fX}{Fitted frequency table}
\item{ResX}{Residual frequency table}
\item{fit}{Fit statistic: weighted least squares
  of deviations between observed and expected frequencies}
\item{item.sc}{Estimated item parameters}
\item{person.sc}{
    Estimated person parameters
            }
\item{ll}{Log-likelihood of the model}
\item{freq.fitted}{Fitted frequencies in a long data frame}
}

\references{
Scheiblechner, H. (1995). Isotonic ordinal
probabilistic models (ISOP). \emph{Psychometrika,
60}, 281-304.

Scheiblechner, H. (1999). Additive conjoint isotonic
probabilistic models (ADISOP). \emph{Psychometrika,
64}, 295-316.
}



%\author{
%Alexander Robitzsch
%}



\note{
For fitting the ADISOP model it is recommended to first
fit the ISOP model and then proceed with the fitted frequency
table from ISOP (see Examples).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For fitting the ISOP model to dichotomous and
polytomous data see \code{\link{isop.dich}}.
}


\examples{
#############################################################################
# EXAMPLE 1: Dataset Reading
#############################################################################

data(data.read)
dat <- as.matrix( data.read)
dat.resp <- 1 - is.na(dat) # response indicator matrix
I <- ncol(dat)

#***
# (1) Data preparation
#     actually only freq.correct and wgt are needed
#     but these matrices must be computed in advance.

# different scores of students
stud.p <- rowMeans( dat, na.rm=TRUE )
# different item p values
item.p <- colMeans( dat, na.rm=TRUE )
item.ps <- sort( item.p, index.return=TRUE)
dat <- dat[,  item.ps$ix ]
# define score groups students
scores <- sort( unique( stud.p ) )
SC <- length(scores)
# create table
freq.correct <- matrix( NA, SC, I )
wgt <- freq.correct
# percent correct
a1 <- stats::aggregate( dat==1, list( stud.p ), mean, na.rm=TRUE )
freq.correct <- a1[,-1]
# weights
a1 <- stats::aggregate( dat.resp, list( stud.p ), sum, na.rm=TRUE )
wgt <- a1[,-1]

#***
# (2) Fit ISOP model
res.isop <- sirt::fit.isop( freq.correct, wgt )
# fitted frequency table
res.isop$fX

#***
# (3) Fit ADISOP model
# use monotonely smoothed frequency table from ISOP model
res.adisop <- sirt::fit.adisop( freq.correct=res.isop$fX, wgt )
# fitted frequency table
res.adisop$fX
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ISOP model}
%% \keyword{ADISOP model}
