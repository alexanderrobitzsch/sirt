%% File Name: xxirt_createThetaDistribution.Rd
%% File Version: 0.11

\name{xxirt_createThetaDistribution}
\alias{xxirt_createThetaDistribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a User Defined Theta Distribution
}

\description{
Creates a user defined theta distribution.
}

\usage{
xxirt_createThetaDistribution(par, est, P, prior=NULL, prior_par1=NULL,
       prior_par2=NULL)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{
Parameter vector with starting values
}
  \item{est}{
Vector of logicals indicating which parameters should be estimated
}
  \item{P}{
Distribution function for \eqn{\bold{\theta}}
}
  \item{prior}{
Prior distribution
}
  \item{prior_par1}{
First parameter of prior distribution
}
  \item{prior_par2}{
Second parameter of prior distribution
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Alexander Robitzsch
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{xxirt}}
}


\examples{
#############################################################################
## EXAMPLE 1: Definition of theta distribution
#############################################################################

#** theta grid
Theta <- matrix( seq(-10,10,length=31), ncol=1 )

#** theta distribution
P_Theta1 <- function( par, Theta, G){
    mu <- par[1]
    sigma <- max( par[2], .01 )
    TP <- nrow(Theta)
    pi_Theta <- matrix( 0, nrow=TP, ncol=G)
    pi1 <- stats::dnorm( Theta[,1], mean=mu, sd=sigma )
    pi1 <- pi1 / sum(pi1)
    pi_Theta[,1] <- pi1
    return(pi_Theta)
                }
#** create distribution class
par_Theta <- c( "mu"=0, "sigma"=1 )
customTheta  <- sirt::xxirt_createThetaDistribution( par=par_Theta,
                       est=c(FALSE,TRUE), P=P_Theta1 )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
