%% File Name: lc.2raters.Rd
%% File Version: 0.17

\name{lc.2raters}
\alias{lc.2raters}
\alias{summary.lc.2raters}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Latent Class Model for Two Exchangeable Raters and One Item
}


\description{
This function computes a latent class model for ratings on an item
based on exchangeable raters (Uebersax & Grove, 1990). Additionally,
several measures of rater agreement are computed (see e.g. Gwet, 2010).
}

\usage{
lc.2raters(data, conv=0.001, maxiter=1000, progress=TRUE)

\method{summary}{lc.2raters}(object,...)
}



%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data frame with item responses (must be ordered from 0 to \eqn{K}) and two
columns which correspond to ratings of two (exchangeable) raters.
}
  \item{conv}{
Convergence criterion
}
  \item{maxiter}{
Maximum number of iterations
}
  \item{progress}{
An optional logical indicating whether iteration progress should be displayed.
}
\item{object}{Object of class \code{lc.2raters}}
\item{\dots}{Further arguments to be passed}
}


\details{
For two exchangeable raters which provide ratings on an item, a latent
class model with \eqn{K+1} classes (if there are \eqn{K+1} item categories
\eqn{0,...,K}) is defined. Where \eqn{P(X=x, Y=y | c)} denotes
the probability that the first rating is \eqn{x} and the second rating is
\eqn{y} given the true but unknown item category (class) \eqn{c}. Ratings are
assumed to be locally independent, i.e.
\deqn{ P(X=x, Y=y | c )=P( X=x | c) \cdot P(Y=y | c )=p_{x|c} \cdot p_{y|c}}
Note that \eqn{P(X=x|c)=P(Y=x|c)=p_{x|c}} holds due to the exchangeability of raters.
The latent class model estimates true class proportions \eqn{\pi_c} and
conditional item probabilities \eqn{p_{x|c}}.
}

\value{
A list with following entries

\item{classprob.1rater.like}{Classification probability \eqn{P(c|x)} of latent
category \eqn{c} given a manifest rating \eqn{x} (estimated by maximum likelihood)}
\item{classprob.1rater.post}{Classification probability \eqn{P(c|x)} of latent
category \eqn{c} given a manifest rating \eqn{x} (estimated by the posterior
distribution)}
\item{classprob.2rater.like}{Classification probability \eqn{P(c|(x,y))}
of latent category \eqn{c} given two manifest ratings \eqn{x} and \eqn{y}
(estimated by maximum likelihood)}
\item{classprob.2rater.post}{Classification probability \eqn{P(c|(x,y))}
of latent category \eqn{c} given two manifest ratings \eqn{x} and \eqn{y}
(estimated by posterior distribution)}
\item{f.yi.qk}{Likelihood of each pair of ratings}
\item{f.qk.yi}{Posterior of each pair of ratings}
\item{probs}{Item response probabilities \eqn{p_{x|c}}}
\item{pi.k}{Estimated class proportions \eqn{\pi_c}}
\item{pi.k.obs}{Observed manifest class proportions}
\item{freq.long}{Frequency table of ratings in long format}
\item{freq.table}{Symmetrized frequency table of ratings}
\item{agree.stats}{Measures of rater agreement. These measures include
percentage agreement (\code{agree0}, \code{agree1}), Cohen's kappa and weighted
Cohen's kappa (\code{kappa}, \code{wtd.kappa.linear}),
Gwet's AC1 agreement measures (\code{AC1}; Gwet, 2008, 2010) and
Aickin's alpha (\code{alpha.aickin}; Aickin, 1990).
}
\item{data}{Used dataset}
\item{N.categ}{Number of categories}
}

\references{
Aickin, M. (1990). Maximum likelihood estimation of agreement in the constant
predictive probability model, and its relation to Cohen's kappa.
\emph{Biometrics}, \bold{46}, 293-302.

Gwet, K. L. (2008). Computing inter-rater reliability and its variance
in the presence of high agreement.
\emph{British Journal of Mathematical and Statistical Psychology},
\bold{61}, 29-48.

Gwet, K. L. (2010). \emph{Handbook of Inter-Rater Reliability}.
Advanced Analytics, Gaithersburg. \url{http://www.agreestat.com/}

Uebersax, J. S., & Grove, W. M. (1990). Latent class analysis of diagnostic
agreement. \emph{Statistics in Medicine}, \bold{9}, 559-572.
}


\author{
Alexander Robitzsch
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{rm.facets}} and \code{\link{rm.sdt}} for
specifying rater models.

See also the \pkg{irr} package for measures of rater agreement.
}
% rm.hrm, rm.facets


\examples{
#############################################################################
# EXAMPLE 1: Latent class models for rating datasets data.si05
#############################################################################

data(data.si05)

#*** Model 1: one item with two categories
mod1 <- sirt::lc.2raters( data.si05$Ex1)
summary(mod1)

#*** Model 2: one item with five categories
mod2 <- sirt::lc.2raters( data.si05$Ex2)
summary(mod2)

#*** Model 3: one item with eight categories
mod3 <- sirt::lc.2raters( data.si05$Ex3)
summary(mod3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{summary}
\keyword{Rater model}
\keyword{Latent class model}
