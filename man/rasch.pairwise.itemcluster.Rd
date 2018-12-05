%% File Name: rasch.pairwise.itemcluster.Rd
%% File Version: 0.37

\name{rasch.pairwise.itemcluster}
\alias{rasch.pairwise.itemcluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pairwise Estimation of the Rasch Model for Locally Dependent Items
}
\description{
This function uses pairwise conditional likelihood estimation
for estimating item parameters in the Rasch model.
}
\usage{
rasch.pairwise.itemcluster(dat, itemcluster=NULL, b.fixed=NULL,
    conv=1e-05, maxiter=3000, progress=TRUE, b.init=NULL, zerosum=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{An \eqn{N \times I} data frame.
  Missing responses are allowed and must be recoded as \code{NA}.
}
  \item{itemcluster}{
Optional integer vector of itemcluster (see Examples). Different
integers correspond to different item clusters. No item cluster
is set as default.
}
  \item{b.fixed}{
Matrix for fixing item parameters. The first columns contains the item
(number or name), the second column the parameter to be fixed.
}
  \item{conv}{
Convergence criterion in maximal absolute parameter change
}
  \item{maxiter}{
Maximal number of iterations
}
  \item{progress}{
A logical which displays progress. Default is \code{TRUE}.
}
  \item{b.init}{
Vector of initial item difficulty estimates. Default is \code{NULL}.
}
\item{zerosum}{Optional logical indicating whether item difficulties
    should be centered in each iteration. The default is that
        no centering is conducted.}

}
\details{
This is an adaptation of the algorithm of van der Linden and Eggen (1986). Only item pairs
of different item clusters are taken into account for item difficulty estimation.
Therefore, the problem of locally dependent items within each itemcluster is (almost)
eliminated (see Examples below) because contributions of local dependencies
do not appear in the pairwise likelihood terms. In detail, the estimation rests
on observed frequency tables of items \eqn{i} and \eqn{j} and therefore on conditional
probabilities
\deqn{ \frac{P(X_i=x, X_j=y)}{P(X_i + X_j=1 )} \quad \mbox{with}
\quad x,y=0,1 \quad \mbox{and} \quad x+y=1 }
If for some item pair \eqn{(i,j)} a higher positive (or negative) correlation
is expected (i.e. deviation from local dependence), then this pair is
removed from estimation. Clearly, there is a loss in precision but item
parameters can be less biased.
}
\value{
Object of class \code{rasch.pairwise} with elements
  \item{b}{
    Vector of item difficulties
         }
 \item{item}{
        Data frame of item parameters (\eqn{N}, \eqn{p} and item difficulty)
            }
}
\references{
van der Linden, W. J., & Eggen, T. J. H. M. (1986). \emph{An empirical
Bayes approach to item banking}. Research Report 86-6, University of Twente.

Zwinderman, A. H. (1995). Pairwise parameter estimation in Rasch models.
\emph{Applied Psychological Measurement, 19}, 369-375.
}



%\author{
%Alexander Robitzsch
%}


\note{
No standard errors are provided by this function. Use resampling
methods for conducting statistical inference.

Formulas for asymptotic standard errors of this pairwise estimation method are
described in Zwinderman (1995).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{rasch.pairwise}}, \code{\link{summary.rasch.pairwise}},

Pairwise marginal likelihood estimation (also labeled as pseudolikelihood
estimation) can be conducted with \code{\link{rasch.pml3}}.

Other estimation methods are implemented in \code{\link{rasch.copula2}} or
\code{\link{rasch.mml2}}.

For simulation of locally dependent data see \code{\link{sim.rasch.dep}}.
}


\examples{
#############################################################################
# EXAMPLE 1: Example with locally dependent items
#      12 Items: Cluster 1 -> Items 1,...,4
#                Cluster 2 -> Items 6,...,9
#                Cluster 3 -> Items 10,11,12
#############################################################################

set.seed(7896)
I <- 12                             # number of items
n <- 5000                           # number of persons
b <- seq(-2,2, len=I)               # item difficulties
bsamp <- b <- sample(b)             # sample item difficulties
theta <- stats::rnorm( n, sd=1 ) # person abilities
# itemcluster
itemcluster <- rep(0,I)
itemcluster[ 1:4 ] <- 1
itemcluster[ 6:9 ] <- 2
itemcluster[ 10:12 ] <- 3
# residual correlations
rho <- c( .55, .25, .45 )

# simulate data
dat <- sirt::sim.rasch.dep( theta, b, itemcluster, rho )
colnames(dat) <- paste("I", seq(1,ncol(dat)), sep="")

# estimation with pairwise Rasch model
mod3 <- sirt::rasch.pairwise( dat )
summary(mod3)

# use item cluster in rasch pairwise estimation
mod <- sirt::rasch.pairwise.itemcluster( dat=dat, itemcluster=itemcluster )
summary(mod)

\dontrun{
# Rasch MML estimation
mod4 <- sirt::rasch.mml2( dat )
summary(mod4)

# Rasch Copula estimation
mod5 <- sirt::rasch.copula2( dat, itemcluster=itemcluster )
summary(mod5)

# compare different item parameter estimates
M1 <- cbind( "true.b"=bsamp,  "b.rasch"=mod4$item$b,
        "b.rasch.copula"=mod5$item$thresh,
        "b.rasch.pairwise"=mod3$b,
        "b.rasch.pairwise.cluster"=mod$b
      )
# center item difficulties
M1 <- scale( M1, scale=FALSE )
round( M1, 3 )
round( apply( M1, 2, stats::sd ), 3 )

#  Below the output of the example is presented.
#  It is surprising that the rasch.pairwise.itemcluster is pretty close
#  to the estimate in the Rasch copula model.

  ##   > M1 <- scale( M1, scale=F )
  ##   > round( M1, 3 )
  ##       true.b b.rasch b.rasch.copula b.rasch.pairwise b.rasch.pairwise.cluster
  ##   I1   0.545   0.561          0.526            0.628                    0.524
  ##   I2  -0.182  -0.168         -0.174           -0.121                   -0.156
  ##   I3  -0.909  -0.957         -0.867           -0.971                   -0.899
  ##   I4  -1.636  -1.726         -1.625           -1.765                   -1.611
  ##   I5   1.636   1.751          1.648            1.694                    1.649
  ##   I6   0.909   0.892          0.836            0.898                    0.827
  ##   I7  -2.000  -2.134         -2.020           -2.051                   -2.000
  ##   I8  -1.273  -1.355         -1.252           -1.303                   -1.271
  ##   I9  -0.545  -0.637         -0.589           -0.581                   -0.598
  ##   I10  1.273   1.378          1.252            1.308                    1.276
  ##   I11  0.182   0.241          0.226            0.109                    0.232
  ##   I12  2.000   2.155          2.039            2.154                    2.026
  ##   > round( apply( M1, 2, sd ), 3 )
  ##                     true.b                  b.rasch           b.rasch.copula
  ##                      1.311                    1.398                    1.310
  ##      b.rasch.pairwise    b.rasch.pairwise.cluster
  ##                 1.373                       1.310   }

# set item parameters of first item to 0 and of second item to -0.7
b.fixed <- cbind( c(1,2), c(0,-.7) )
mod5 <- sirt::rasch.pairwise.itemcluster( dat=dat, b.fixed=b.fixed,
             itemcluster=itemcluster )
# difference between estimations 'mod' and 'mod5'
dfr <- cbind( mod5$item$b, mod$item$b )
plot( mod5$item$b, mod$item$b, pch=16)
apply( dfr, 1, diff )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Pairwise conditional maximum likelihood (PCML)}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\keyword{Local dependence}
