%% File Name: btm.Rd
%% File Version: 0.31

\name{btm}
\alias{btm}
\alias{summary.btm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extended Bradley-Terry Model
}
\description{
Estimates an extended Bradley-Terry model (Hunter, 2004; see Details).
Parameter estimation uses a bias corrected joint maximum likelihood
estimation method based on
\eqn{\varepsilon}-adjustment (see Bertoli-Barsotti, Lando & Punzo, 2014).
See Details for the algorithm.
}

\usage{
btm(data, ignore.ties=FALSE, fix.eta=NULL, fix.delta=NULL, fix.theta=NULL,
       maxiter=100, conv=1e-04, eps=0.3)

\method{summary}{btm}(object, file=NULL, digits=4,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data frame with three columns. The first two columns contain labels
from the units in the pair comparison. The third column contains the
result of the comparison. "1" means that the first units wins, "0" means
that the second unit wins and "0.5" means a draw (a tie).
}
  \item{ignore.ties}{
Logical indicating whether ties should be ignored.
}
  \item{fix.eta}{
Numeric value for a fixed \eqn{\eta} value
}
  \item{fix.delta}{
Numeric value for a fixed \eqn{\delta} value
}
  \item{fix.theta}{
A vector with entries for fixed theta values.
}
  \item{maxiter}{
Maximum number of iterations
}
  \item{conv}{
Convergence criterion
}
  \item{eps}{
The \eqn{\varepsilon} parameter for the \eqn{\varepsilon}-adjustment
method (see Bertoli-Barsotti, Lando & Punzo, 2014) which reduces bias
in ability estimates. In case of \eqn{\varepsilon=0}, persons with
extreme scores are removed from the pairwise comparison.
}
\item{object}{Object of class \code{btm}}
\item{file}{Optional file name for sinking the summary into}
\item{digits}{Number of digits after decimal to print}
\item{\dots}{Further arguments to be passed.}
}

\details{
The extended Bradley-Terry model for the comparison of individuals
\eqn{i} and \eqn{j} is defined as
\deqn{P(X_{ij}=1 ) \propto \exp( \eta + \theta_i ) }
\deqn{P(X_{ij}=0 ) \propto \exp(  \theta_j ) }
\deqn{P(X_{ij}=0.5) \propto \exp( \delta + ( \eta  + \theta_i +\theta_j )/2 ) }
The parameters \eqn{\theta_i} denote the abilities, \eqn{\delta} is the
tendency of the occurrence of ties and \eqn{\eta} is the home-advantage
effect.

A joint maximum likelihood (JML) estimation is applied for simulataneous 
estimation of \eqn{\eta}, \eqn{\delta} and all \eqn{\theta_i} parameters.
In the Rasch model, it was shown that JML can result in biased parameter
estimates. The \eqn{\varepsilon}-adjustment approach has been proposed 
to reduce the bias in parameter estimates (Bertoli-Bersotti, Lando & Punzo, 2014).
This estimation approach is adapted to the Bradley-Terry model in
the \code{btm} function. To this end, the likelihood function is
modified for the purpose of bias reduction. It can be easily shown that there 
exist sufficient statistics for \eqn{\eta}, \eqn{\delta} and all \eqn{\theta_i} 
parameters. In the \eqn{\varepsilon}-adjustment approach, the sufficient
statistic for the \eqn{\theta_i} parameter is modified. In JML estimation
of the Bradley-Terry model, \eqn{S_i= \sum_{j \ne i} ( x_{ij} + x_{ji} )} is
a suifficient statistic for \eqn{\theta_i}. Let \eqn{M_i} the maximum score
for person \eqn{i} which is the number of \eqn{x_{ij}} terms appaering in
\eqn{S_i}. In the \eqn{\varepsilon}-adjustment approach, the sufficient statistic
\eqn{S_i} is modified to 
\deqn{S_{i, \varepsilon} = \varepsilon + 
\frac{M_i - 2 \varepsilon}{M_i} S_i } and \eqn{S_{i, \varepsilon}} intead of
\eqn{S_{i}} is used in JML estimation. Hence, original scores \eqn{S_i} are
linearly transformed for all persons \eqn{i}.
}

\value{
List with following entries
  \item{pars}{Parameter summary for \eqn{\eta} and \eqn{\delta}}
  \item{effects}{Parameter estimates for \eqn{\theta} and
      outfit and infit statistics}
  \item{summary.effects}{Summary of \eqn{\theta} parameter estimates}
  \item{mle.rel}{MLE reliability, also known as separation reliability}
  \item{sepG}{Separation index \eqn{G}}
  \item{probs}{Estimated probabilities}
  \item{data}{Used dataset with integer identifiers}
}

\references{
Bertoli-Barsotti, L., Lando, T., & Punzo, A. (2014). Estimating a Rasch Model 
via fuzzy empirical probability functions. In D. Vicari, A. Okada, G. Ragozini 
& C. Weihs (Eds.). \emph{Analysis and Modeling of Complex Data in Behavioral 
and Social Sciences}. Springer. 

Hunter, D. R. (2004). MM algorithms for generalized Bradley-Terry models.
\emph{Annals of Statistics, 32}, 384-406.
}

\author{
Alexander Robitzsch
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also the \R packages \pkg{BradleyTerry2}, \pkg{psychotools},
\pkg{psychomix} and \pkg{prefmod}.
}


\examples{
#############################################################################
# EXAMPLE 1: Bradley-Terry model | data.pw01
#############################################################################

data(data.pw01)

dat <- data.pw01
dat <- dat[, c("home_team", "away_team", "result") ]

# recode results according to needed input
dat$result[ dat$result==0 ] <- 1/2   # code for ties
dat$result[ dat$result==2 ] <- 0     # code for victory of away team

#********************
# Model 1: Estimation with ties and home advantage
mod1 <- sirt::btm( dat)
summary(mod1)

\dontrun{
#********************
# Model 2: Estimation with ties, no epsilon adjustment
mod2 <- sirt::btm( dat, eps=0, fix.eta=0)
summary(mod2)

#********************
# Model 3: Some fixed abilities
fix.theta <- c("Anhalt Dessau"=-1 )
mod3 <- sirt::btm( dat, eps=0, fix.theta=fix.theta)
summary(mod3)

#********************
# Model 4: Ignoring ties, no home advantage effect
mod4 <- sirt::btm( dat, ignore.ties=TRUE, fix.eta=0)
summary(mod4)

#********************
# Model 5: Ignoring ties, no home advantage effect (JML approach -> eps=0)
mod5 <- sirt::btm( dat, ignore.ties=TRUE, fix.eta=0, eps=0)
summary(mod5)

#############################################################################
# EXAMPLE 2: Venice chess data
#############################################################################

# See http://www.rasch.org/rmt/rmt113o.htm
# Linacre, J. M. (1997). Paired Comparisons with Standard Rasch Software.
# Rasch Measurement Transactions, 11:3, 584-585.

# dataset with chess games -> "D" denotes a draw (tie)
chessdata <- scan( what="character")
    1D.0..1...1....1.....1......D.......D........1.........1.......... Browne
    0.1.D..0...1....1.....1......D.......1........D.........1......... Mariotti
    .D0..0..1...D....D.....1......1.......1........1.........D........ Tatai
    ...1D1...D...D....1.....D......D.......D........1.........0....... Hort
    ......010D....D....D.....1......D.......1........1.........D...... Kavalek
    ..........00DDD.....D.....D......D.......1........D.........1..... Damjanovic
    ...............00D0DD......D......1.......1........1.........0.... Gligoric
    .....................000D0DD.......D.......1........D.........1... Radulov
    ............................DD0DDD0D........0........0.........1.. Bobotsov
    ....................................D00D00001.........1.........1. Cosulich
    .............................................0D000D0D10..........1 Westerinen
    .......................................................00D1D010000 Zichichi

L <- length(chessdata) / 2
games <- matrix( chessdata, nrow=L, ncol=2, byrow=TRUE )
G <- nchar(games[1,1])
# create matrix with results
results <- matrix( NA, nrow=G, ncol=3 )
for (gg in 1:G){
    games.gg <- substring( games[,1], gg, gg )
    ind.gg <- which( games.gg !="." )
    results[gg, 1:2 ] <- games[ ind.gg, 2]
    results[gg, 3 ] <- games.gg[ ind.gg[1] ]
            }
results <- as.data.frame(results)
results[,3] <- paste(results[,3] )
results[ results[,3]=="D", 3] <- 1/2
results[,3] <- as.numeric( results[,3] )

# fit model ignoring draws
mod1 <- sirt::btm( results, ignore.ties=TRUE, fix.eta=0, eps=0 )
summary(mod1)

# fit model with draws
mod2 <- sirt::btm( results, fix.eta=0, eps=0 )
summary(mod2)
}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Bradley-Terry model}
\keyword{summary}% __ONLY ONE__ keyword per line
