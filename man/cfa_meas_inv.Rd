%% File Name: cfa_meas_inv.Rd
%% File Version: 0.074

\name{cfa_meas_inv}
\alias{cfa_meas_inv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of a Unidimensional Factor Model under Full and Partial
Measurement Invariance
}
\description{
Estimates a unidimensional factor model based on the normal distribution fitting
function under full and partial measurement invariance.
Item loadings and item intercepts are successively freed based on the largest
modification index and a chosen significance level \code{alpha}.
}
\usage{
cfa_meas_inv(dat, group, weights=NULL, alpha=0.01, verbose=FALSE, op=c("~1","=~"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame containing items
}
  \item{group}{
Vector of group identifiers
}
  \item{weights}{
Optional vector of sampling weights
}
  \item{alpha}{
Significance level
}
  \item{verbose}{
Logical indicating whether progress should be shown
}
\item{op}{Operators (intercepts or loadings) for which estimation should be freed}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
List with several entries
  \item{pars_mi}{Model parameters under full invariance}
  \item{pars_pi}{Model parameters under partial invariance}
  \item{mod_mi}{Fitted model under full invariance}
  \item{mod_pi}{Fitted model under partial invariance}

  \item{\ldots}{More output}
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link[sirt:invariance.alignment]{sirt::invariance.alignment}}
}
\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Factor model under full and partial invariance
#############################################################################

#--- data simulation

set.seed(65)
G <- 3  # number of groups
I <- 5  # number of items
# define lambda and nu parameters
lambda <- matrix(1, nrow=G, ncol=I)
nu <- matrix(0, nrow=G, ncol=I)
err_var <- matrix(1, nrow=G, ncol=I)

# define size of noninvariance
dif <- 1
#- 1st group: N(0,1)
lambda[1,3] <- 1+dif*.4; nu[1,5] <- dif*.5
#- 2nd group: N(0.3,1.5)
gg <- 2 ;
lambda[gg,5] <- 1-.5*dif; nu[gg,1] <- -.5*dif
#- 3nd group: N(.8,1.2)
gg <- 3
lambda[gg,4] <- 1-.7*dif; nu[gg,2] <- -.5*dif
#- define distributions of groups
mu <- c(0,.3,.8)
sigma <- sqrt(c(1,1.5,1.2))
N <- rep(1000,3) # sample sizes per group

#* use simulation function
dat <- sirt::invariance_alignment_simulate(nu, lambda, err_var, mu, sigma, N,
                exact=TRUE)

#--- estimate CFA
mod <- sirt::cfa_meas_inv(dat=dat[,-1], group=dat$group, verbose=TRUE, alpha=0.05)
mod$pars_mi
mod$pars_pi
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
