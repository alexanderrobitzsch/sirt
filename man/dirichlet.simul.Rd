%% File Name: dirichlet.simul.Rd
%% File Version: 0.10

\name{dirichlet.simul}
\alias{dirichlet.simul}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of a Dirichlet Distributed Vectors
}
\description{
This function makes random draws from a Dirichlet distribution.
}
\usage{
dirichlet.simul(alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{
A matrix with \eqn{\bold{\alpha}} parameters of the
Dirichlet distribution
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A data frame with Dirichlet distributed responses
}
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

%\seealso{
%For examples see \code{\link{dirichlet.mle}}
%}

\examples{
#############################################################################
# EXAMPLE 1: Simulation with two components
#############################################################################

set.seed(789)
N <- 2000
probs <- c(.7 , .3)	# define (extremal) class probabilities

#*** alpha0 = .2  -> nearly crisp latent classes
alpha0 <- .2
alpha <- alpha0*probs
alpha <- matrix( alpha , nrow=N , ncol=length(alpha) , byrow=TRUE  )
x <- sirt::dirichlet.simul( alpha )
htitle <- expression(paste( alpha[0], "=.2, ", p[1] , "=.7"   ) )
hist( x[,1] , breaks = seq(0,1,len=20) , main=htitle)

#*** alpha0 = 3 -> strong deviation from crisp membership
alpha0 <- 3
alpha <- alpha0*probs
alpha <- matrix( alpha , nrow=N , ncol=length(alpha) , byrow=TRUE  )
x <- sirt::dirichlet.simul( alpha )
htitle <- expression(paste( alpha[0], "=3, ", p[1] , "=.7"   ) )
hist( x[,1] , breaks = seq(0,1,len=20) , main=htitle)

\dontrun{
#############################################################################
# EXAMPLE 2: Simulation with three components
#############################################################################

set.seed(986)
N <- 2000
probs <- c( .5 , .35 , .15 )

#*** alpha0 = .2
alpha0 <- .2
alpha <- alpha0*probs
alpha <- matrix( alpha , nrow=N , ncol=length(alpha) , byrow=TRUE  )
x <- sirt::dirichlet.simul( alpha )
htitle <- expression(paste( alpha[0], "=.2, ", p[1] , "=.7"   ) )
miceadds::library_install("ade4")
ade4::triangle.plot(x, label=NULL , clabel = 1)

#*** alpha0 = 3
alpha0 <- 3
alpha <- alpha0*probs
alpha <- matrix( alpha , nrow=N , ncol=length(alpha) , byrow=TRUE  )
x <- sirt::dirichlet.simul( alpha )
htitle <- expression(paste( alpha[0], "=3, ", p[1] , "=.7"   ) )
ade4::triangle.plot(x, label=NULL , clabel = 1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Dirichlet distribution}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
