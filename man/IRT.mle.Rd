\name{IRT.mle}
\alias{IRT.mle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Person Parameter Estimation
}
\description{
Computes the maximum likelihood estimate (MLE), 
weighted likelihood estimate (WLE) and maximum aposterior
estimate (MAP) of ability in unidimensional item response models 
(Penfield & Bergeron, 2005; Warm, 1989). Item response functions can be 
defined by the user.
}

\usage{
IRT.mle(data, irffct, arg.list, theta=rep(0,nrow(data)), type = "MLE", 
     mu=0, sigma=1, maxiter = 20, maxincr = 3, h = 0.001, convP = 1e-04, 
     maxval = 9, progress = TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data frame with item responses
}
\item{irffct}{
User defined item response (see Examples). Arguments must be
specified in \code{arg.list}. The function must contain \code{theta}
and \code{ii} (item index) as arguments.
}
  \item{theta}{
Initial ability estimate
}
  \item{arg.list}{
List of arguments for \code{irffct}.
}
  \item{type}{Type of ability estimate. It can be \code{"MLE"} (the default),
  	  \code{"WLE"} or \code{"MAP"}. 
}
\item{mu}{Mean of normal prior distriubution (for \code{type="MAP"})}
\item{sigma}{Standard deviation of normal prior distriubution (for \code{type="MAP"})}
  \item{maxiter}{
Maximum number of iterations
}
  \item{maxincr}{
Maximum increment
}
  \item{h}{
Numerical differentiation parameter
}
  \item{convP}{
Convergence criterion
}
  \item{maxval}{
Maximum ability value to be estimated
}
  \item{progress}{
Logical indicating whether iteration progress should be displayed
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
Data frame with estimated abilities (\code{est}) and its standard error
(\code{se}).
}

\references{
Penfield, R. D., & Bergeron, J. M. (2005). Applying a weighted 
maximum likelihood latent trait estimator to the generalized 
partial credit model. \emph{Applied Psychological Measurement},
\bold{29}, 218-233.	
	
Warm, T. A. (1989). Weighted likelihood estimation of ability in item response theory.
\emph{Psychometrika}, \bold{54}, 427-450.
}

\author{
Alexander Robitzsch
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also the \pkg{PP} package for further person parameter
estimation methods.
}

\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Generalized partial credit model
#############################################################################

data(data.ratings1)
dat <- data.ratings1

# estimate model
mod1 <- sirt::rm.facets( dat[ , paste0( "k",1:5) ], rater=dat$rater, 
             pid=dat$idstud  , maxiter=15)
# extract dataset and item parameters
data <- mod1$procdata$dat2.NA
a <- mod1$ipars.dat2$a
b <- mod1$ipars.dat2$b
theta0 <- mod1$person$EAP
# define item response function for item ii
calc.pcm <- function( theta , a , b , ii ){
    K <- base::ncol(b)
    N <- base::length(theta)
    matrK <- base::matrix( 0:K , nrow=N , ncol=K+1 , byrow=TRUE)
    eta <- a[ii] * theta * matrK - base::matrix( c(0,b[ii,]), nrow=N, ncol=K+1, byrow=TRUE)
    eta <- base::exp(eta)
    probs <- eta / base::rowSums(eta, na.rm=TRUE)    
    base::return(probs)
}
arg.list <- list("a"=a , "b"=b )

# MLE
abil1 <- sirt::IRT.mle( data, irffct=calc.pcm, theta=theta0, arg.list=arg.list )
str(abil1)
# WLE
abil2 <- sirt::IRT.mle( data, irffct=calc.pcm, theta=theta0, arg.list=arg.list, type="WLE")
str(abil2)
# MAP with prior distribution N(.2, 1.3)
abil3 <- sirt::IRT.mle( data, irffct=calc.pcm, theta=theta0, arg.list=arg.list, 
              type="MAP", mu=.2, sigma=1.3 )
str(abil3)

#############################################################################
# EXAMPLE 2: Rasch model
#############################################################################

data(data.read)
dat <- data.read
I <- ncol(dat)

# estimate Rasch model
mod1 <- sirt::rasch.mml2( dat )
summary(mod1)

# define item response function
irffct <- function( theta, b , ii){
    eta <- exp( theta - b[ii] )
    probs <- eta / ( 1 + eta )
    probs <- cbind( 1 - probs , probs )
    return(probs)
}
# initial person parameters and item parameters
theta0 <- mod1$person$EAP
arg.list <- list( "b" = mod1$item$b  )

# estimate WLE
abil <- sirt::IRT.mle( data = dat , irffct=irffct , arg.list=arg.list , 
            theta=theta0, type="WLE")
# compare with wle.rasch function
theta <- sirt::wle.rasch( dat , b= mod1$item$b )
cbind( abil[,1] , theta$theta , abil[,2] , theta$se.theta )

#############################################################################
# EXAMPLE 3: Ramsay quotient model
#############################################################################

data(data.read)
dat <- data.read
I <- ncol(dat)

# estimate Ramsay model
mod1 <- sirt::rasch.mml2( dat , irtmodel ="ramsay.qm" )
summary(mod1)
# define item response function
irffct <- function( theta, b , K , ii){
    eta <- exp( theta / b[ii] )
    probs <- eta / ( K[ii] + eta )
    probs <- cbind( 1 - probs , probs )
    return(probs)
}
# initial person parameters and item parameters
theta0 <- exp( mod1$person$EAP )
arg.list <- list( "b" = mod1$item2$b , "K"=mod1$item2$K )
# estimate MLE
res <- sirt::IRT.mle( data = dat , irffct=irffct , arg.list=arg.list , theta=theta0 ,
            maxval=20 , maxiter=50)
}	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{Person parameters}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
