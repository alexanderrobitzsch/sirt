%% File Name: mle.pcm.group.Rd
%% File Version: 0.22

\name{mle.pcm.group}
\alias{mle.pcm.group}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum Likelihood Estimation of Person or Group Parameters
in the Generalized Partial Credit Model
}
\description{
This function estimates person or group parameters in the
partial credit model (see Details).
}
\usage{
mle.pcm.group(dat, b, a=rep(1, ncol(dat)), group=NULL,
    pid=NULL, adj_eps=0.3, conv=1e-04, maxiter=30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
A numeric \eqn{N \times I} matrix
}
  \item{b}{
Matrix with item thresholds
}
  \item{a}{
Vector of item slopes
}
  \item{group}{
Vector of group identifiers
}
  \item{pid}{
Vector of person identifiers
}
  \item{adj_eps}{
Numeric value which is used in \eqn{\varepsilon} adjustment
of the likelihood. A value of zero (or a very small
\eqn{\varepsilon>0}) corresponds to the usual maximum
likelihood estimate.
}
  \item{conv}{
Convergence criterion
}
  \item{maxiter}{
Maximum number of iterations
}
}
\details{
It is assumed that the generalized partial credit model holds.
In case one estimates a person parameter \eqn{\theta_p},
the log-likelihood is maximized and the following
estimating equation results: (see Penfield & Bergeron, 2005):
\deqn{ 0=( \log L )'=\sum_i a_i \cdot [ \tilde{x}_{pi} -
    E(X_{pi} | \theta_p ) ] }
where \eqn{E(X_{pi} | \theta_p )} denotes the expected item response
conditionally on \eqn{\theta_p}.

With the method of \eqn{\varepsilon}-adjustment (Bertoli-Barsotti & Punzo, 2012;
Bertoli-Barsotti, Lando & Punzo, 2014),
the observed item responses \eqn{x_{pi}} are transformed such that
no perfect scores arise and bias is reduced. If \eqn{S_p}
is the sum score of person \eqn{p} and \eqn{M_p} the maximum
score of this person, then the transformed sum scores \eqn{\tilde{S}_p}
are
\deqn{ \tilde{S}_p=\varepsilon + \frac{M_p - 2 \varepsilon}{M_p} S_p}
However, the adjustment is directly conducted on item responses
to also handle the case of the generalized partial credit
model with item slope parameters different from 1.

In case one estimates a group parameter \eqn{\theta_g},
the following estimating equation is used:
\deqn{ 0=(\log L )'=\sum_p \sum_i a_i \cdot [ \tilde{x}_{pgi} -
    E(X_{pgi} | \theta_g ) ] }
where persons \eqn{p} are nested within a group \eqn{g}.
The \eqn{\varepsilon}-adjustment is then performed at the
group level, not at the individual level.
}
\value{
A list with following entries:
  \item{person}{Data frame with person or group parameters}
  \item{data_adjeps}{Modified dataset according to the
  \eqn{\varepsilon} adjustment.
  }
}
\references{
Bertoli-Barsotti, L., & Punzo, A. (2012).
Comparison of two bias reduction techniques for the Rasch model.
\emph{Electronic Journal of Applied Statistical Analysis, 5}, 360-366.

Bertoli-Barsotti, L., Lando, T., & Punzo, A. (2014).
Estimating a Rasch Model via fuzzy empirical probability functions.
In D. Vicari, A. Okada, G. Ragozini &  C. Weihs (Eds.).
\emph{Analysis and Modeling of Complex Data in Behavioral and Social Sciences},
Springer.

Penfield, R. D., & Bergeron, J. M. (2005). Applying a weighted
maximum likelihood latent trait estimator to the generalized
partial credit model. \emph{Applied Psychological Measurement,
29}, 218-233.
}
\author{
Alexander Robitzsch
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Estimation of a group parameter for only one item per group
#############################################################################

data(data.si01)
dat <- data.si01
# item parameter estimation (partial credit model) in TAM
library(TAM)
mod <- TAM::tam.mml( dat[,2:3], irtmodel="PCM")
# extract item difficulties
b <- matrix( mod$xsi$xsi, nrow=2, byrow=TRUE )
# groupwise estimation
res1 <- sirt::mle.pcm.group( dat[,2:3], b=b, group=dat$idgroup )
# individual estimation
res2 <- sirt::mle.pcm.group( dat[,2:3], b=b  )

#############################################################################
# EXAMPLE 2: Data Reading data.read
#############################################################################

data(data.read)
# estimate Rasch model
mod <- sirt::rasch.mml2( data.read )
score <- rowSums( data.read )
data.read <- data.read[ order(score), ]
score <- score[ order(score) ]
# compare different epsilon-adjustments
res30 <- sirt::mle.pcm.group( data.read, b=matrix( mod$item$b, 12, 1 ),
               adj_eps=.3 )$person
res10 <- sirt::mle.pcm.group( data.read, b=matrix( mod$item$b, 12, 1 ),
             adj_eps=.1 )$person
res05 <- sirt::mle.pcm.group( data.read, b=matrix( mod$item$b, 12, 1 ),
              adj_eps=.05 )$person
# plot different scorings
plot( score, res05$theta, type="l", xlab="Raw score", ylab=expression(theta[epsilon]),
         main="Scoring with different epsilon-adjustments")
lines( score, res10$theta, col=2, lty=2 )
lines( score, res30$theta, col=4, lty=3 )
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Person parameters}
\keyword{Group parameters}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
