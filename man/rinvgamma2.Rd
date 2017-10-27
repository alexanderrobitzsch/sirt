%% File Name: rinvgamma2.Rd
%% File Version: 0.11

\name{rinvgamma2}
\alias{rinvgamma2}
\alias{dinvgamma2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inverse Gamma Distribution in Prior Sample Size Parameterization
}
\description{
Random draws and density of inverse gamma distribution parametrized
in prior sample size \code{n0} and prior variance \code{var0}
(see Gelman et al., 2014).
}
\usage{
rinvgamma2(n, n0, var0)

dinvgamma2(x, n0, var0)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of draws for inverse gamma distribution
}
  \item{n0}{
Prior sample size
}
  \item{var0}{
Prior variance
}
\item{x}{Vector with numeric values for density evaluation}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A vector containing random draws or density values	
}


%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...


\references{
Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., 
& Rubin, D. B. (2014). 
\emph{Bayesian data analysis} (Vol. 3). Boca Raton, FL, USA: Chapman & Hall/CRC.
}

\author{
Alexander Robitzsch
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{MCMCpack::rinvgamma},
\code{\link[stats:GammaDist]{stats::rgamma}},
\code{MCMCpack::dinvgamma},
\code{\link[stats:GammaDist]{stats::dgamma}}
}
\examples{
#############################################################################
# EXAMPLE 1: Inverse gamma distribution 
#############################################################################

# prior sample size of 100 and prior variance of 1.5
n0 <- 100
var0 <- 1.5

# 100 random draws
y1 <- sirt::rinvgamma2( n=100 , n0 , var0 )
summary(y1)
graphics::hist(y1)

# density y at grid x
x <- seq( 0 , 2 , len=100 )
y <- sirt::dinvgamma2( x , n0 , var0 )
graphics::plot( x , y , type="l")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
