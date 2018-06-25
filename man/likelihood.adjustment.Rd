%% File Name: likelihood.adjustment.Rd
%% File Version: 0.10

\name{likelihood.adjustment}
\alias{likelihood.adjustment}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Adjustment and Approximation of Individual Likelihood Functions
}

\description{
Approximates individual likelihood functions \eqn{L(\bold{X}_p | \theta)}
by normal distributions (see Mislevy, 1990). Extreme response patterns
are handled by adding pseudo-observations of items with extreme item
difficulties (see argument \code{extreme.item}. The individual standard
deviations of the likelihood, used in the normal approximation, can be
modified by individual adjustment factors which are specified in \code{adjfac}.
In addition, a reliability of the adjusted likelihood can be specified
in \code{target.EAP.rel}.
}

\usage{
likelihood.adjustment(likelihood, theta=NULL, prob.theta=NULL,
     adjfac=rep(1, nrow(likelihood)), extreme.item=5, target.EAP.rel=NULL,
     min_tuning=0.2, max_tuning=3, maxiter=100, conv=1e-04,
     trait.normal=TRUE)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{likelihood}{
A matrix containing the individual likelihood \eqn{L(\bold{X}_p | \theta)} or
an object of class \code{IRT.likelihood}.
}
  \item{theta}{
Optional vector of (unidimensional) \eqn{\theta} values
}
  \item{prob.theta}{
Optional vector of probabilities of \eqn{\theta} trait distribution
}
  \item{adjfac}{
Vector with individual adjustment factors of the standard deviations of the
likelihood
}
  \item{extreme.item}{
Item difficulties of two extreme pseudo items which are added as additional
observed data to the likelihood. A large number (e.g. \code{extreme.item=15})
leaves the likelihood almost unaffected. See also Mislevy (1990).
}
  \item{target.EAP.rel}{
Target EAP reliability. An additional tuning parameter is estimated
which adjusts the likelihood to obtain a pre-specified reliability.
}
  \item{min_tuning}{
Minimum value of tuning parameter (if \code{! is.null(target.EAP.rel) })
}
  \item{max_tuning}{
Maximum value of tuning parameter (if \code{! is.null(target.EAP.rel) })
}
  \item{maxiter}{
Maximum number of iterations (if \code{! is.null(target.EAP.rel) })
}
  \item{conv}{
Convergence criterion (if \code{! is.null(target.EAP.rel) })
}
  \item{trait.normal}{
Optional logical indicating whether the trait distribution should be
normally distributed (if \code{! is.null(target.EAP.rel) }).
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
Object of class \code{IRT.likelihood}.
}

\references{
Mislevy, R. (1990). Scaling procedures. In E. Johnson & R. Zwick (Eds.),
\emph{Focusing the new design: The NAEP 1988 technical report} (ETS RR 19-20).
Princeton, NJ: Educational Testing Service.
}


\author{
Alexander Robitzsch
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[CDM:IRT.likelihood]{CDM::IRT.likelihood}},
\code{\link[TAM:tam.latreg]{TAM::tam.latreg}}
}

\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Adjustment of the likelihood | data.read
#############################################################################

library(CDM)
library(TAM)
data(data.read)
dat <- data.read

# define theta grid
theta.k <- seq(-6,6,len=41)

#*** Model 1: fit Rasch model in TAM
mod1 <- TAM::tam.mml( dat, control=list( nodes=theta.k) )
summary(mod1)

#*** Model 2: fit Rasch copula model
testlets <- substring( colnames(dat), 1, 1 )
mod2 <- sirt::rasch.copula2( dat, itemcluster=testlets, theta.k=theta.k)
summary(mod2)

# model comparison
IRT.compareModels( mod1, mod2 )

# extract EAP reliabilities
rel1 <- mod1$EAP.rel
rel2 <- mod2$EAP.Rel
# variance inflation factor
vif <- (1-rel2) / (1-rel1)
  ##  > vif
  ##  [1] 1.211644

# extract individual likelihood
like1 <- IRT.likelihood( mod1 )
# adjust likelihood from Model 1 to obtain a target EAP reliability of .599
like1b <- sirt::likelihood.adjustment( like1, target.EAP.rel=.599 )

# compare estimated latent regressions
lmod1a <- TAM::tam.latreg( like1, Y=NULL )
lmod1b <- TAM::tam.latreg( like1b, Y=NULL )
summary(lmod1a)
summary(lmod1b)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Likelihood adjustment}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
