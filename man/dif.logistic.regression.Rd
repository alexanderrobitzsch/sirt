%% File Name: dif.logistic.regression.Rd
%% File Version: 0.24

\name{dif.logistic.regression}
\alias{dif.logistic.regression}
\title{
Differential Item Functioning using Logistic Regression Analysis
}
\description{
This function estimates differential item
functioning using a logistic regression analysis
(Zumbo, 1999). 
}
\usage{
dif.logistic.regression(dat, group, score,quant=1.645)
}
\arguments{
  \item{dat}{
Data frame with dichotomous item responses
}
  \item{group}{
Group identifier
}
  \item{score}{
Ability estimate, e.g. the WLE.
}
\item{quant}{Used quantile of the normal distribution 
for assessing statistical significance}
}

\details{
Items are classified into A (negligible DIF), B (moderate DIF) and 
C (large DIF) levels according to the
ETS classification system (Longford, Holland & Thayer, 1993, p. 175).
See also Monahan et al. (2007) for further DIF effect size
classifications.
}

\value{
A data frame with following variables:

\item{itemnr}{Numeric index of the item}
\item{sortDIFindex}{Rank of item with respect to the uniform DIF
	(from negative to positive values)
}
\item{item}{Item name}
\item{N}{Sample size per item}
\item{R}{Value of \code{group} variable for reference group}
\item{F}{Value of \code{group} variable for focal group}
\item{nR}{Sample size per item in reference group}
\item{nF}{Sample size per item in focal group}
\item{p}{Item \eqn{p} value}
\item{pR}{Item \eqn{p} value in reference group}
\item{pF}{Item \eqn{p} value in focal group}
\item{pdiff}{Item \eqn{p} value differences}
\item{pdiff.adj}{Adjusted \eqn{p} value difference}
\item{uniformDIF}{Uniform DIF estimate}
\item{se.uniformDIF}{Standard error of uniform DIF}
\item{t.uniformDIF}{The \eqn{t} value for uniform DIF}
\item{sig.uniformDIF}{Significance label for uniform DIF}
\item{DIF.ETS}{DIF classification according to the ETS classification
	system (see Details)
}
\item{uniform.EBDIF}{Empirical Bayes estimate of uniform DIF (Longford,
	Holland & Thayer, 1993) which takes degree of DIF standard error
        into account}
\item{DIF.SD}{Value of the DIF standard deviation}
\item{nonuniformDIF}{Nonuniform DIF estimate}
\item{se.nonuniformDIF}{Standard error of nonuniform DIF}
\item{t.nonuniformDIF}{The \eqn{t} value for nonuniform DIF}
\item{sig.nonuniformDIF}{Significance label for nonuniform DIF}
}


\references{
Longford, N. T., Holland, P. W., & Thayer, D. T. (1993).
Stability of the MH D-DIF statistics across populations.
In P. W. Holland & H. Wainer (Eds.). \emph{Differential
Item Functioning} (pp. 171-196). Hillsdale, NJ: Erlbaum.	

Monahan, P. O., McHorney, C. A., Stump, T. E., & Perkins, A. J. (2007). 
Odds ratio, delta, ETS classification, and standardization measures of 
DIF magnitude for binary logistic regression. 
\emph{Journal of Educational and Behavioral Statistics}, \bold{32}, 92-109.

Zumbo, B. D. (1999). \emph{A handbook on the theory and methods of differential 
item functioning (DIF): Logistic regression modeling as a unitary framework 
for binary and Likert-type (ordinal) item scores}. 
Ottawa ON: Directorate of Human Resources Research and Evaluation, 
Department of National Defense. 
}
\author{
Alexander Robitzsch
}
\seealso{
For assessing DIF variance see \code{\link{dif.variance}} and
\code{\link{dif.strata.variance}}

See also \code{\link{rasch.evm.pcm}} for assessing differential item
functioning in the partial credit model.

See the \pkg{difR} package for a large collection of DIF detection
methods.

For a download of the free \emph{DIF-Pack} software (SIBTEST, ...) see  
\url{http://psychometrictools.measuredprogress.org/home}.

}
\examples{
#############################################################################
# EXAMPLE 1: Mathematics data | Gender DIF
#############################################################################

data( data.math )
dat <- data.math$data
items <- grep( "M" , colnames(dat))

# estimate item parameters and WLEs
mod <- sirt::rasch.mml2( dat[,items] )
wle <- sirt::wle.rasch( dat[,items] , b=mod$item$b )$theta

# assess DIF by logistic regression
mod1 <- sirt::dif.logistic.regression( dat=dat[,items] , score=wle , group=dat$female)

# calculate DIF variance
dif1 <- sirt::dif.variance( dif=mod1$uniformDIF , se.dif = mod1$se.uniformDIF )
dif1$unweighted.DIFSD
  ## > dif1$unweighted.DIFSD
  ## [1] 0.1963958

# calculate stratified DIF variance
# stratification based on domains
dif2 <- sirt::dif.strata.variance( dif=mod1$uniformDIF , se.dif = mod1$se.uniformDIF ,
              itemcluster = data.math$item$domain )
  ## $unweighted.DIFSD
  ## [1] 0.1455916

\dontrun{
#****
# Likelihood ratio test and graphical model test in eRm package
miceadds::library_install("eRm")
# estimate Rasch model
res <- eRm::RM( dat[,items] )
summary(res)
# LR-test with respect to female
lrres <- eRm::LRtest(res, splitcr = dat$female)
summary(lrres)
# graphical model test
eRm::plotGOF(lrres) 

#############################################################################
# EXAMPLE 2: Comparison with Mantel-Haenszel test
#############################################################################

library(TAM)
library(difR)

#*** (1) simulate data
set.seed(776)
N <- 1500   # number of persons per group
I <- 12     # number of items
mu2 <- .5   # impact (group difference)
sd2 <- 1.3  # standard deviation group 2

# define item difficulties
b <- seq( -1.5 , 1.5 , length=I)
# simulate DIF effects
bdif <- scale( stats::rnorm(I , sd = .6 ) , scale=FALSE )[,1]
# item difficulties per group
b1 <- b + 1/2 * bdif
b2 <- b - 1/2 * bdif
# simulate item responses
dat1 <- sirt::sim.raschtype( theta = stats::rnorm(N , mean=0 , sd =1 ) , b = b1 )
dat2 <- sirt::sim.raschtype( theta = stats::rnorm(N , mean=mu2 , sd = sd2 ) , b = b2 )
dat <- rbind( dat1 , dat2 )
group <- rep( c(1,2) , each=N ) # define group indicator

#*** (2) scale data 
mod <- TAM::tam.mml( dat , group=group )
summary(mod)

#*** (3) extract person parameter estimates
mod_eap <- mod$person$EAP
mod_wle <- tam.wle( mod )$theta

#*********************************
# (4) techniques for assessing differential item functioning

# Model 1: assess DIF by logistic regression and WLEs
dif1 <- sirt::dif.logistic.regression( dat=dat , score= mod_wle , group= group)
# Model 2: assess DIF by logistic regression and EAPs
dif2 <- sirt::dif.logistic.regression( dat=dat , score= mod_eap , group= group)
# Model 3: assess DIF by Mantel-Haenszel statistic
dif3 <- difR::difMH(Data=dat, group=group, focal.name="1" ,  purify=FALSE )
print(dif3)
  ##  Mantel-Haenszel Chi-square statistic: 
  ##   
  ##        Stat.    P-value     
  ##  I0001  14.5655   0.0001 ***
  ##  I0002 300.3225   0.0000 ***
  ##  I0003   2.7160   0.0993 .  
  ##  I0004 191.6925   0.0000 ***
  ##  I0005   0.0011   0.9740    
  ##  [...]
  ##  Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1  
  ##  Detection threshold: 3.8415 (significance level: 0.05)
  ##
  ##  Effect size (ETS Delta scale): 
  ##   
  ##  Effect size code: 
  ##   'A': negligible effect 
  ##   'B': moderate effect 
  ##   'C': large effect 
  ##   
  ##        alphaMH deltaMH  
  ##  I0001  1.3908 -0.7752 A
  ##  I0002  0.2339  3.4147 C
  ##  I0003  1.1407 -0.3093 A
  ##  I0004  2.8515 -2.4625 C
  ##  I0005  1.0050 -0.0118 A
  ##  [...]
  ##  
  ##  Effect size codes: 0 'A' 1.0 'B' 1.5 'C' 
  ##   (for absolute values of 'deltaMH') 

# recompute DIF parameter from alphaMH
uniformDIF3 <- log(dif3$alphaMH) 

# compare different DIF statistics
dfr <- data.frame( "bdif"=  bdif , "LR_wle" = dif1$uniformDIF ,
        "LR_eap" = dif2$uniformDIF , "MH" = uniformDIF3 )
round( dfr , 3 )
  ##       bdif LR_wle LR_eap     MH
  ##  1   0.236  0.319  0.278  0.330
  ##  2  -1.149 -1.473 -1.523 -1.453
  ##  3   0.140  0.122  0.038  0.132
  ##  4   0.957  1.048  0.938  1.048
  ##  [...]
colMeans( abs( dfr[,-1] - bdif ))
  ##      LR_wle     LR_eap         MH 
  ##  0.07759187 0.19085743 0.07501708 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Differential item functioning (DIF)}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
