%% File Name: Q3.Rd
%% File Version: 1.13

\name{Q3}
\alias{Q3}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of the \eqn{Q_3} Statistic (Yen, 1984)
}
\description{
This function estimates the \eqn{Q_3} statistic according to Yen (1984).
The statistic \eqn{Q_3} is calculated for every item pair \eqn{(i,j)}
which is the correlation between item residuals after fitting the Rasch model.
}
\usage{
Q3(dat, theta, b, progress=TRUE)
}
\arguments{
  \item{dat}{
An \eqn{N \times I} data frame of dichotomous item responses
}
  \item{theta}{
Vector of length \eqn{N} of person parameter estimates (e.g. obtained from
    \code{\link{wle.rasch}})
}
  \item{b}{
Vector of length \eqn{I} (e.g. obtained from \code{\link{rasch.mml2}})
}
  \item{progress}{
Should iteration progress be displayed?
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with following entries
\item{q3.matrix}{An \eqn{I \times I} matrix of \eqn{Q_3} statistics}
\item{q3.long}{Just the \code{q3.matrix} in long matrix format where every row
                  corresponds to an item pair}
\item{expected}{An \eqn{N \times I} matrix of expected probabilities
    by the Rasch model}
\item{residual}{An \eqn{N \times I} matrix of residuals obtained after
    fitting the Rasch model}
\item{Q3.stat}{Vector with descriptive statistics of \eqn{Q_3}}
%% ...
%%        res <- list( "q3.matrix" = q3.matr , "q3.long" = q3.long ,
%%        "expected" = expected , "residual" = residual )
}
\references{
Yen, W. M. (1984). Effects of local item dependence on the fit and equating
performance of the three-parameter logistic model.
\emph{Applied Psychological Measurement}, \bold{8}, 125-145.
}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For the estimation of the average \eqn{Q_3} statistic within testlets see
\code{\link{Q3.testlet}}.

For modeling testlet effects see \code{\link{mcmc.3pno.testlet}}.

For handling local dependencies in IRT models see
\code{\link{rasch.copula2}}, \code{\link{rasch.pml3}} or \cr
\code{\link{rasch.pairwise.itemcluster}}.
}
\examples{
#############################################################################
# EXAMPLE 1: data.read. The 12 items are arranged in 4 testlets
#############################################################################
data(data.read)

# estimate the Rasch model
mod <- sirt::rasch.mml2( data.read)
# estmate WLEs
mod.wle <- sirt::wle.rasch( dat = data.read , b = mod$item$b )
# calculate Yen's Q3 statistic
mod.q3 <- sirt::Q3( dat = data.read , theta = mod.wle$theta , b = mod$item$b )
  ##   Yen's Q3 Statistic based on an estimated theta score
  ##   *** 12 Items | 66 item pairs
  ##   *** Q3 Descriptives
  ##        M     SD    Min    10%    25%    50%    75%    90%    Max
  ##   -0.085  0.110 -0.261 -0.194 -0.152 -0.107 -0.051  0.041  0.412

# plot Q3 statistics
I <- ncol(data.read)
image( 1:I , 1:I , mod.q3$q3.matrix , col=gray( 1 - (0:32)/32) ,
        xlab="Item" , ylab="Item")
abline(v=c(5,9)) # borders for testlets
abline(h=c(5,9))

\dontrun{
# obtain Q3 statistic from modelfit.sirt function which is based on the
# posterior distribution of theta and not on observed values
fitmod <- sirt::modelfit.sirt( mod )
# extract Q3 statistic
q3stat <- fitmod$itempairs$Q3
  ##  > summary(q3stat)
  ##      Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
  ##  -0.21760 -0.11590 -0.07280 -0.05545 -0.01220  0.44710
  ##  > sd(q3stat)
  ##  [1] 0.1101451
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Testlets}
\keyword{Local dependence}
\keyword{Q3}

