%% File Name: pgenlogis.Rd
%% File Version: 0.27

\name{pgenlogis}
\alias{pgenlogis}
\alias{genlogis.moments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of Probabilities and Moments for the
Generalized Logistic Item Response Model
}
\description{
Calculation of probabilities and moments for the generalized logistic
item response model (Stukel, 1988).
}
\usage{
pgenlogis(x, alpha1=0, alpha2=0)

genlogis.moments(alpha1, alpha2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vector
}
  \item{alpha1}{
Upper tail parameter \eqn{\alpha_1} in the generalized
logistic item response model. The default is 0.
}
  \item{alpha2}{
Lower tail parameter \eqn{\alpha_2} parameter in the generalized
logistic item response model. The default is 0.
}
}
\details{
The class of generalized logistic link functions contain
the most important link functions using the specifications (Stukel, 1988):

\itemize{
\item logistic link function \eqn{L}:
\deqn{ L(x) \approx G_{ ( \alpha_1=0, \alpha_2=0)}[ x ] }

\item probit link function \eqn{\Phi}:
\deqn{ \Phi(x) \approx G_{ ( \alpha_1=0.165, \alpha_2=0.165)}[ 1.47 x ] }


\item loglog link function \eqn{H}:
\deqn{ H(x) \approx G_{ (\alpha_1=-0.037, \alpha_2=0.62)}[
-0.39+1.20x-0.007x^2] }

\item cloglog link function \eqn{H}:
\deqn{ H(x) \approx G_{ ( \alpha_1=0.62, \alpha_2=-0.037)}[
0.54+1.64x+0.28x^2+0.046x^3] }

}
}

% # loglog link
% y <- pgenlogis( -.39 + 1.20*x -.007*x^2, alpha1=-.037, alpha2=.62 )

% # cloglog link
% y <- pgenlogis( .54+1.64*x +.28*x^2 + .046*x^3, alpha1=.062, alpha2=-.037 )

\value{
Vector of probabilities or moments
}

\references{
Stukel, T. A. (1988). Generalized logistic models.
\emph{Journal of the American Statistical Association,
83}, 426-431.
}


%\author{
%Alexander Robitzsch
%}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
sirt::pgenlogis( x=c(-.3, 0, .25, 1 ), alpha1=0, alpha2=.6 )
  ##   [1] 0.4185580 0.5000000 0.5621765 0.7310586

####################################################################
# compare link functions
x <- seq( -3,3, .1 )

#***
# logistic link
y <- sirt::pgenlogis( x, alpha1=0, alpha2=0 )
plot( x, stats::plogis(x), type="l", main="Logistic Link", lwd=2)
points( x, y, pch=1, col=2 )

#***
# probit link
round( sirt::genlogis.moments( alpha1=.165, alpha2=.165 ), 3 )
  ##       M    SD   Var
  ##   0.000 1.472 2.167
# SD of generalized logistic link function is 1.472
y <- sirt::pgenlogis( x * 1.47, alpha1=.165, alpha2=.165 )
plot( x, stats::pnorm(x), type="l", main="Probit Link", lwd=2)
points( x, y, pch=1, col=2 )

#***
# loglog link
y <- sirt::pgenlogis( -.39 + 1.20*x -.007*x^2, alpha1=-.037, alpha2=.62 )
plot( x, exp( - exp( -x ) ), type="l", main="Loglog Link", lwd=2,
    ylab="loglog(x)=exp(-exp(-x))" )
points( x, y, pch=17, col=2 )

#***
# cloglog link
y <- sirt::pgenlogis( .54+1.64*x +.28*x^2 + .046*x^3, alpha1=.062, alpha2=-.037 )
plot( x, 1-exp( - exp(x) ), type="l", main="Cloglog Link", lwd=2,
    ylab="loglog(x)=1-exp(-exp(x))" )
points( x, y, pch=17, col=2 )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Generalized logistic item response model}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
