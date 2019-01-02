%% File Name: wle.rasch.Rd
%% File Version: 1.18

\name{wle.rasch}
\alias{wle.rasch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weighted Likelihood Estimation of Person Abilities
}
\description{
This function computes weighted likelihood estimates for dichotomous responses based
on the Rasch model (Warm, 1989).
}
\usage{
wle.rasch(dat, dat.resp=NULL, b, itemweights=1 + 0 * b,
    theta=rep(0, nrow(dat)), conv=0.001, maxit=200,
    wle.adj=0, progress=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{ An \eqn{N \times I} data frame of
  dichotomous item responses
}
  \item{dat.resp}{
    Optional data frame with dichotomous response indicators
}
  \item{b}{
    Vector of length \eqn{I} with fixed item difficulties
}
  \item{itemweights}{
    Optional vector of fixed item discriminations
}
  \item{theta}{
    Optional vector of initial person parameter estimates
}
  \item{conv}{
    Convergence criterion
}
  \item{maxit}{
    Maximal number of iterations
}
  \item{wle.adj}{
    Constant for WLE adjustment
}
  \item{progress}{
    Display progress?
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
A list with following entries
  \item{theta}{
        Estimated weighted likelihood estimate
        }
  \item{dat.resp}{
        Data frame with dichotomous response indicators. A one indicates
        an observed response, a zero a missing response. See also \code{dat.resp}
        in the list of arguments of this function.
        }
   \item{p.ia}{Matrix with expected item response, i.e.
    the probabilities \eqn{P(X_{pi}=1|\theta_p )=invlogit( \theta_p - b_i )}.
            }
   \item{wle}{WLE reliability (Adams, 2005)}
}


\references{
Adams, R. J. (2005). Reliability as a measurement design effect.
\emph{Studies in Educational Evaluation, 31}, 162-172.

Warm, T. A. (1989). Weighted likelihood estimation of ability in item response theory.
\emph{Psychometrika, 54}, 427-450.
}



%\author{
%Alexander Robitzsch
%}



%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For standard errors of weighted likelihood estimates estimated via jackknife
see \code{\link{wle.rasch.jackknife}}.

For a joint estimation of item and person parameters see the joint maximum
likelihood estimation method in \code{\link{rasch.jml}}.
}
\examples{
#############################################################################
# EXAMPLE 1: Dataset Reading
#############################################################################
data(data.read)

# estimate the Rasch model
mod <- sirt::rasch.mml2(data.read)
mod$item

# estmate WLEs
mod.wle <- sirt::wle.rasch( dat=data.read, b=mod$item$b )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Weighted likelihood estimation (WLE)}
%% \keyword{Person parameter estimation}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
