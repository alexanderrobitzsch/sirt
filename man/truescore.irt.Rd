%% File Name: truescore.irt.Rd
%% File Version: 0.15

\name{truescore.irt}
\alias{truescore.irt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conversion of Trait Scores \eqn{\theta} into
True Scores \eqn{\tau ( \theta )}
}
\description{
This function computes the true score
\eqn{\tau=\tau(\theta) = \sum_{i=1}^I P_i(\theta)}
in a unidimensional item response model with \eqn{I} items. In addition, it also
transforms conditional standard errors if they are provided.
}
\usage{
truescore.irt(A, B, c = NULL, d = NULL, theta = seq(-3, 3, len = 21),
    error = NULL, pid = NULL, h = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
Matrix or vector of item slopes. See Examples for polytomous responses.
}
  \item{B}{
Matrix or vector of item intercepts. Note that the entries in
\code{B} refer to item intercepts and not to item difficulties.
}
  \item{c}{
Optional vector of guessing parameters
}
  \item{d}{
Optional vector of slipping parameters
}
  \item{theta}{
Vector of trait values
}
  \item{error}{
Optional vector of standard errors of trait
}
  \item{pid}{
Optional vector of person identifiers
}
  \item{h}{
Numerical differentiation parameter
}
}

\details{
In addition, the function \eqn{\pi(\theta)= \frac{1}{I} \cdot \tau( \theta)}
of the expected percent score is approximated by a logistic function
\deqn{ \pi ( \theta ) \approx l + ( u - l ) \cdot invlogit ( a \theta + b )
}
}

\value{
A data frame with following columns:
\item{truescore}{True scores \eqn{\tau = \tau ( \theta )}}
\item{truescore.error}{Standard errors of true scores}
\item{percscore}{Expected correct scores which is \eqn{\tau}
divided by the maximum true score}
\item{percscore.error}{Standard errors of expected correct scores}
\item{lower}{The \eqn{l} parameter}
\item{upper}{The \eqn{u} parameter}
\item{a}{The \eqn{a} parameter}
\item{b}{The \eqn{b} parameter}
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

%}

\examples{
#############################################################################
# EXAMPLE 1: Dataset with mixed dichotomous and polytomous responses
#############################################################################

data(data.mixed1)
dat <- data.mixed1

#****
# Model 1: Partial credit model
# estimate model with TAM package
library(TAM)
mod1 <- TAM::tam.mml( dat )
# estimate person parameter estimates
wmod1 <- TAM::tam.wle( mod1 )
wmod1 <- wmod1[ order(wmod1$theta) , ]
# extract item parameters
A <- mod1$B[,-1,1]
B <- mod1$AXsi[,-1]
# person parameters and standard errors
theta <- wmod1$theta
error <- wmod1$error

# estimate true score transformation
dfr <- sirt::truescore.irt( A=A , B=B , theta=theta , error=error )

# plot different person parameter estimates and standard errors
par(mfrow=c(2,2))
plot( theta , dfr$truescore , pch=16 , cex=.6 , xlab=expression(theta) , type="l",
    ylab=expression(paste( tau , "(",theta , ")" ))  , main="True Score Transformation" )
plot( theta , dfr$percscore , pch=16 , cex=.6 , xlab=expression(theta) , type="l",
    ylab=expression(paste( pi , "(",theta , ")" ))  , main="Percent Score Transformation" )
points( theta , dfr$lower + (dfr$upper-dfr$lower)*
                stats::plogis(dfr$a*theta+dfr$b) , col=2 , lty=2)
plot( theta , error , pch=16 , cex=.6 , xlab=expression(theta) , type="l",
    ylab=expression(paste("SE(",theta , ")" ))  , main="Standard Error Theta" )
plot( dfr$truescore , dfr$truescore.error , pch=16 , cex=.6 , xlab=expression(tau) ,
    ylab=expression(paste("SE(",tau , ")" ) ) , main="Standard Error True Score Tau" ,
    type="l")
par(mfrow=c(1,1))

\dontrun{
#****
# Model 2: Generalized partial credit model
mod2 <- TAM::tam.mml.2pl( dat , irtmodel="GPCM")
# estimate person parameter estimates
wmod2 <- TAM::tam.wle( mod2 )
# extract item parameters
A <- mod2$B[,-1,1]
B <- mod2$AXsi[,-1]
# person parameters and standard errors
theta <- wmod2$theta
error <- wmod2$error
# estimate true score transformation
dfr <- sirt::truescore.irt( A=A , B=B , theta=theta , error=error )

#############################################################################
# EXAMPLE 2: Dataset Reading data.read
#############################################################################
data(data.read)

#****
# Model 1: estimate difficulty + guessing model
mod1 <- sirt::rasch.mml2( data.read , fixed.c = rep(.25,12) )
mod1$person <- mod1$person[ order( mod1$person$EAP) , ]
# person parameters and standard errors
theta <- mod1$person$EAP
error <- mod1$person$SE.EAP
A <- rep(1,12)
B <- - mod1$item$b
c <- rep(.25,12)
# estimate true score transformation
dfr <- sirt::truescore.irt( A=A , B=B , theta=theta , error=error ,c=c)

plot( theta , dfr$percscore , pch=16 , cex=.6 , xlab=expression(theta) , type="l",
    ylab=expression(paste( pi , "(",theta , ")" ))  , main="Percent Score Transformation" )
points( theta , dfr$lower + (dfr$upper-dfr$lower)*
             stats::plogis(dfr$a*theta+dfr$b) , col=2 , lty=2)

#****
# Model 2: Rasch model
mod2 <- sirt::rasch.mml2( data.read  )
# person parameters and standard errors
theta <- mod2$person$EAP
error <- mod2$person$SE.EAP
A <- rep(1,12)
B <- - mod2$item$b
# estimate true score transformation
dfr <- sirt::truescore.irt( A=A , B=B , theta=theta , error=error )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{True scores}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
