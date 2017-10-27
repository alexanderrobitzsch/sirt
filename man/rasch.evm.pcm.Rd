%% File Name: rasch.evm.pcm.Rd
%% File Version: 0.33

\name{rasch.evm.pcm}
\alias{rasch.evm.pcm}
\alias{summary.rasch.evm.pcm}
\alias{coef.rasch.evm.pcm}
\alias{vcov.rasch.evm.pcm}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of the Partial Credit Model using the Eigenvector Method
}

\description{
This function performs the eigenvector approach to estimate item
parameters which is based on a pairwise estimation approach
(Gardner & Engelhard, 2002). No assumption about person parameters
is required for item parameter estimation. Statistical inference is performed by 
Jackknifing. If a group identifier is provided, tests for differential item
functioning are performed. 
}

\usage{
rasch.evm.pcm(dat, jackunits = 20, weights = NULL, pid = NULL ,
    group=NULL , powB = 2, adj_eps = 0.3, progress = TRUE )

\method{summary}{rasch.evm.pcm}(object,...)

\method{coef}{rasch.evm.pcm}(object,...)

\method{vcov}{rasch.evm.pcm}(object,...)	
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame with dichotomous or polytomous item responses 
}
\item{jackunits}{
A number of Jackknife units (if an integer is provided as the argument
value) or a vector in which the Jackknife units are already
defined.
}
  \item{weights}{
Optional vector of sample weights
}
  \item{pid}{
Optional vector of person identifiers
}
\item{group}{
Optional vector of group identifiers. In this case, item parameters are group wise
estimated and tests for differential item functioning are performed.
}
  \item{powB}{
Power created in \eqn{B} matrix which is the basis of parameter estimation
}
  \item{adj_eps}{
Adjustment parameter for person parameter estimation
(see \code{\link{mle.pcm.group}})
}
  \item{progress}{
An optional logical indicating whether progress should be displayed
}
\item{object}{Object of class \code{rasch.evm.pcm}}
\item{\dots}{Further arguments to be passed}
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}


\value{
A list with following entries

\item{item}{Data frame with item parameters. The item parameter estimate
	is denoted by \code{est} while a Jackknife bias-corrected estimate
	is \code{est_jack}. The Jackknife standard error is \code{se}.
}
\item{b}{Item threshold parameters}
\item{person}{Data frame with person parameters obtained (MLE)}
\item{B}{Paired comparison matrix}
\item{D}{Transformed paired comparison matrix}
\item{coef}{Vector of estimated coefficients}
\item{vcov}{Covariance matrix of estimated item parameters}
\item{JJ}{Number of jackknife units}
\item{JJadj}{Reduced number of jackknife units}
\item{powB}{Used power of comparison matrix \eqn{B}}
\item{maxK}{Maximum number of categories per item}
\item{G}{Number of groups}
\item{desc}{Some descriptives}
\item{difstats}{Statistics for differential item functioning if \code{group}
is provided as an argument}

}

\references{
Choppin, B. (1985). A fully conditional estimation procedure for Rasch Model
parameters. \emph{Evaluation in Education}, \bold{9}, 29-42.

Garner, M., & Engelhard, G. J. (2002). An eigenvector method for 
estimating item parameters of the dichotomous and polytomous Rasch models. 
\emph{Journal of Applied Measurement}, \bold{3}, 107-128.

Wang, J., & Engelhard, G. (2014). A pairwise algorithm in \R for rater-mediated
assessments. \emph{Rasch Measurement Transactions}, \bold{28(1)}, 1457-1459.
}

\author{
Alexander Robitzsch
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See the \pkg{pairwise} package for the alternative row averaging
approach of Choppin (1985) and Wang and Engelhard (2014) for an
alternative \R implementation.
}


\examples{
#############################################################################
# EXAMPLE 1: Dataset Liking for Science
#############################################################################

data(data.liking.science)
dat <- data.liking.science

# estimate partial credit model using 10 Jackknife units
mod1 <- sirt::rasch.evm.pcm( dat , jackunits=10 )
summary(mod1)

\dontrun{
# compare results with TAM
library(TAM)
mod2 <- TAM::tam.mml( dat )
r1 <- mod2$xsi$xsi
r1 <- r1 - mean(r1)
# item parameters are similar
dfr <- data.frame( "b_TAM"=r1 , mod1$item[,c( "est","est_jack") ] )
round( dfr , 3 )
  ##      b_TAM    est est_jack
  ##  1  -2.496 -2.599   -2.511
  ##  2   0.687  0.824    1.030
  ##  3  -0.871 -0.975   -0.943
  ##  4  -0.360 -0.320   -0.131
  ##  5  -0.833 -0.970   -0.856
  ##  6   1.298  1.617    1.444
  ##  7   0.476  0.465    0.646
  ##  8   2.808  3.194    3.439
  ##  9   1.611  1.460    1.433
  ##  10  2.396  1.230    1.095
  ##  [...]

# partial credit model in eRm package
miceadds::library_install("eRm")
mod3 <- eRm::PCM(X=dat)
summary(mod3)
eRm::plotINFO(mod3)      # plot item and test information
eRm::plotICC(mod3)       # plot ICCs
eRm::plotPImap(mod3)     # plot person-item maps

#############################################################################
# EXAMPLE 2: Garner and Engelhard (2002) toy example dichotomous data
#############################################################################

dat <- scan()
   1 0 1 1   1 1 0 0   1 0 0 0   0 1 1 1   1 1 1 0   
   1 1 0 1   1 1 1 1   1 0 1 0   1 1 1 1   1 1 0 0

dat <- matrix( dat , 10 , 4 , byrow=TRUE)
colnames(dat) <- paste0("I" , 1:4 )

# estimate Rasch model with no jackknifing
mod1 <- sirt::rasch.evm.pcm( dat , jackunits=0 )

# paired comparison matrix
mod1$B
  ##          I1_Cat1 I2_Cat1 I3_Cat1 I4_Cat1
  ##  I1_Cat1       0       3       4       5
  ##  I2_Cat1       1       0       3       3
  ##  I3_Cat1       1       2       0       2
  ##  I4_Cat1       1       1       1       0

#############################################################################
# EXAMPLE 3: Garner and Engelhard (2002) toy example polytomous data
#############################################################################

dat <- scan()
   2 2 1 1 1   2 1 2 0 0   1 0 0 0 0   0 1 1 2 0   1 2 2 1 1   
   2 2 0 2 1   2 2 1 1 0   1 0 1 0 0   2 1 2 2 2   2 1 0 0 1

dat <- matrix( dat , 10 , 5 , byrow=TRUE)
colnames(dat) <- paste0("I" , 1:5 )

# estimate partial credit model with no jackknifing
mod1 <- sirt::rasch.evm.pcm( dat , jackunits=0 , powB=3 )

# paired comparison matrix
mod1$B
  ##          I1_Cat1 I1_Cat2 I2_Cat1 I2_Cat2 I3_Cat1 I3_Cat2 I4_Cat1 I4_Cat2 I5_Cat1 I5_Cat2
  ##  I1_Cat1       0       0       2       0       1       1       2       1       2       1
  ##  I1_Cat2       0       0       0       3       2       2       2       2       2       3
  ##  I2_Cat1       1       0       0       0       1       1       2       0       2       1
  ##  I2_Cat2       0       1       0       0       1       2       0       3       1       3
  ##  I3_Cat1       1       1       1       1       0       0       1       2       3       1
  ##  I3_Cat2       0       1       0       2       0       0       1       1       1       1
  ##  I4_Cat1       0       1       0       0       0       2       0       0       1       2
  ##  I4_Cat2       1       0       0       2       1       1       0       0       1       1
  ##  I5_Cat1       0       1       0       1       2       1       1       2       0       0
  ##  I5_Cat2       0       0       0       1       0       0       0       0       0       0

#############################################################################
# EXAMPLE 4: Partial credit model for dataset data.mg from CDM package
#############################################################################

library(CDM)
data(data.mg,package="CDM")
dat <- data.mg[ , paste0("I",1:11) ]

#*** Model 1: estimate partial credit model
mod1 <- sirt::rasch.evm.pcm( dat )
# item parameters
round( mod1$b , 3 )
  ##        Cat1   Cat2  Cat3
  ##  I1  -1.537     NA    NA
  ##  I2  -2.360     NA    NA
  ##  I3  -0.574     NA    NA
  ##  I4  -0.971 -2.086    NA
  ##  I5  -0.104  0.201    NA
  ##  I6   0.470  0.806    NA
  ##  I7  -1.027  0.756 1.969
  ##  I8   0.897     NA    NA
  ##  I9   0.766     NA    NA
  ##  I10  0.069     NA    NA
  ##  I11 -1.122  1.159 2.689

#*** Model 2: estimate PCM with pairwise package
miceadds::library_install("pairwise")
mod2 <- pairwise::pair(daten=dat)
summary(mod2)
plot(mod2)
# compute standard errors
semod2 <- pairwise::pairSE(daten=dat,  nsample = 20)
semod2

#############################################################################
# EXAMPLE 5: Differential item functioning for dataset data.mg
#############################################################################

library(CDM)
data(data.mg,package="CDM")
dat <- data.mg[ data.mg$group \%in\% c(2,3,11) , ]
# define items
items <- paste0("I",1:11)
# estimate model
mod1 <- sirt::rasch.evm.pcm( dat[,items] , weights= dat$weight , group= dat$group )
summary(mod1)

#############################################################################
# EXAMPLE 6: Differential item functioning for Rasch model
#############################################################################

# simulate some data
set.seed(9776)
N <- 1000	# number of persons
I <- 10		# number of items
# simulate data for first group
b <- seq(-1.5,1.5,len=I)
dat1 <- sirt::sim.raschtype( stats::rnorm(N) , b )
# simulate data for second group
b1 <- b
b1[4] <- b1[4] + .5 # introduce DIF for fourth item
dat2 <- sirt::sim.raschtype( stats::rnorm(N,mean=.3) , b1 )
dat <- rbind(dat1 , dat2 )
group <- rep( 1:2 , each=N )
# estimate model
mod1 <- sirt::rasch.evm.pcm( dat , group= group )
summary(mod1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Eigenvector method}
\keyword{Pairwise estimation}
\keyword{summary}
\keyword{coef}
\keyword{vcov}
