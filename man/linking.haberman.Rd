%% File Name: linking.haberman.Rd
%% File Version: 0.755

\name{linking.haberman}
\alias{linking.haberman}
\alias{summary.linking.haberman}
\alias{linking.haberman.lq}
\alias{summary.linking.haberman.lq}
\alias{linking_haberman_itempars_prepare}
\alias{linking_haberman_itempars_convert}
\alias{L0_polish}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linking in the 2PL/Generalized Partial Credit Model
}
\description{
This function does the linking of several studies which are calibrated
using the 2PL or the generalized item response model according to
Haberman (2009). This method is a generalization of log-mean-mean
linking from one study to several studies. The default \code{a_log=TRUE}
logarithmizes item slopes for linking while otherwise an additive regression
model is assumed for the original item loadings (see Details; Battauz, 2017)
}
\usage{
linking.haberman(itempars, personpars, estimation="OLS", a_trim=Inf, b_trim=Inf,
    lts_prop=.5, a_log=TRUE, conv=1e-05, maxiter=1000, progress=TRUE,
    adjust_main_effects=TRUE, vcov=TRUE)

\method{summary}{linking.haberman}(object, digits=3, file=NULL, ...)

linking.haberman.lq(itempars, pow=2, eps=1e-3, a_log=TRUE, use_nu=FALSE,
      est_pow=FALSE, lower_pow=.1, upper_pow=3)

\method{summary}{linking.haberman.lq}(object, digits=3, file=NULL, ...)

## prepare 'itempars' argument for linking.haberman()
linking_haberman_itempars_prepare(b, a=NULL, wgt=NULL)

## conversion of different parameterizations of item parameters
linking_haberman_itempars_convert(itempars=NULL, lambda=NULL, nu=NULL, a=NULL, b=NULL)

## L0 polish precedure minimizing number of interactions in two-way table
L0_polish(x, tol, conv=0.01, maxiter=30, type=1, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{itempars}{
A data frame with four or five columns. The first four columns contain
in the order: study name, item name, \eqn{a} parameter, \eqn{b} parameter.
The fifth column is an optional weight for every item and every study.
}
\item{personpars}{A list with vectors (e.g. EAPs or WLEs) or data frames
    (e.g. plausible values) containing person parameters which
    should be transformed.
If a data frame in each list entry has \code{se} or \code{SE}
(standard error) in a column name, then the corresponding
column is only multiplied by \eqn{A_t}.
If a column is labeled as \code{pid} (person ID),
then it is left untransformed.
}
\item{estimation}{Estimation method. Can be \code{"OLS"} (ordinary least squares),
\code{"BSQ"} (bisquare weighted regression), \code{"HUB"} (regression using Huber
weights), \code{"MED"} (median regression), \code{"LTS"} (trimmed least squares),
\code{"L1"} (median polish), \code{"L0"} (minimizing number of interactions)
}
\item{a_trim}{Trimming parameter for item slopes \eqn{a_{it}} in
bisquare regression (see Details).}
\item{b_trim}{Trimming parameter for item slopes \eqn{b_{it}} in
bisquare regression (see Details).}
\item{lts_prop}{Proportion of retained observations in \code{"LTS"}
regression estimation}
\item{a_log}{Logical indicating whether item slopes should be logarithmized
for linking.}
\item{conv}{
Convergence criterion.
}
  \item{maxiter}{
Maximum number of iterations.
}
  \item{progress}{
An optional logical indicating whether computational progress
should be displayed.
}
\item{adjust_main_effects}{Logical indicating whether all elements in the vector
of main effects should be simultaneously adjusted}
\item{vcov}{Optional indicating whether covariance matrix for linking errors
should be computed}
\item{pow}{Power \eqn{q}}
\item{eps}{Epsilon value used in differentiable approximating function}
\item{use_nu}{Logical indicating whether item intercepts instead of
item difficulties are used in linking}
\item{est_pow}{Logical indicating whether power values should be estimated}
\item{lower_pow}{Lower bound for estimated power}
\item{upper_pow}{Upper bound for estimated power}
\item{lambda}{Matrix containing item loadings}
\item{nu}{Matrix containing item intercepts}
\item{object}{Object of class \code{linking.haberman}.}
\item{digits}{Number of digits after decimals for rounding in \code{summary}.}
\item{file}{Optional file name if \code{summary} should be sunk into a file.}
\item{\dots}{Further arguments to be passed}
\item{b}{Matrix of item intercepts (items \eqn{times} studies)}
\item{a}{Matrix of item slopes}
\item{wgt}{Matrix of weights}
\item{x}{Matrix}
\item{tol}{Tolerance value}
\item{type}{Can be \code{1} (using Tukey's median polish) or
\code{2} (alternating median regression).}
\item{verbose}{Logical indicating whether iteration progress should be displayed}
}

% L0_polish(x, tol, eps=0.01, maxiter=30, type=1, verbose=TRUE)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETAILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



\details{
For \eqn{t=1,\ldots,T} studies, item difficulties \eqn{b_{it}} and
item slopes \eqn{a_{it}} are available. For dichotomous responses, these
parameters are defined by the 2PL response equation
\deqn{ logit P(X_{pi}=1| \theta_p )=a_i ( \theta_p - b_i ) }
while for polytomous responses the generalized partial credit model holds
\deqn{ log \frac{P(X_{pi}=k| \theta_p )}{P(X_{pi}=k-1| \theta_p )}
=a_i ( \theta_p - b_i + d_{ik} ) }

The parameters \eqn{ \{ a_{it}, b_{it} \}} of all items and studies are
linearly transformed using equations \eqn{a_{it} \approx a_i / A_t}
(if \code{a_log=TRUE}) or \eqn{a_{it} \approx a_i + A_t}
(if \code{a_log=FALSE}) and
\eqn{b_{it} \cdot A_t \approx B_t + b_i}. For identification reasons,
we define \eqn{A_1=1} and \eqn{B_1}=0.

The optimization function (which is a least squares criterion;
see Haberman, 2009) seeks the transformation parameters \eqn{A_t} and
\eqn{B_t} with an alternating least squares
method (\code{estimation="OLS"}). Note that every item \eqn{i} and every study \eqn{t} can
be weighted (specified in the fifth column of \code{itempars}).
Alternatively, a robust regression method based on bisquare weighting (Fox, 2015)
can be employed for linking using the argument \code{estimation="BSQ"}.
For example, in the case of item loadings, bisquare weighting is applied to
residuals \eqn{e_{it}=a_{it} - a_i - A_t } (where logarithmized or non-logarithmized
item loadings are employed) forming weights
\eqn{w_{it}=[ 1 - ( e_{it} / k )^2 ]^2} for \eqn{e_{it} <k} and 0 for \eqn{e_{it} \ge k}
where \eqn{k} is the trimming constant which can be estimated or fixed
during estimation using arguments \code{a_trim} or \code{b_trim}. Items in studies with
large residuals
(i.e., presence differential item functioning) are effectively set to zero in the
linking procedure. Alternatively, Huber weights (\code{estimation="HUB"}) downweight
large residuals by applying \eqn{w_{it}=k / | e_{it} |} for residuals
\eqn{|e_{it}|>k}.  The method \code{estimation="LTS"} employs trimmed least squares
where the proportion
of data retained is specified in \code{lts_prop} with default set to .50.

The method \code{estimation="MED"} estimates item parameters and linking constants
based on alternating median regression. A similar approach is the median polish
procedure of Tukey (Tukey, 1977, p. 362ff.; Maronna, Martin & Yohai, 2006, p. 104;
see also \code{\link[stats:medpolish]{stats::medpolish}}) implemented in
\code{estimation="L1"} which aims to minimize \eqn{\sum_{i,t} | e_{it} |}.
For a pre-specified tolerance value \eqn{t} (in \code{a_trim} or \code{b_trim}),
the approach \code{estimation="L0"} minimizes the number of interactions
(i.e., DIF effects) in the \eqn{e_{it}} effects. In more detail, it minimizes
 \eqn{\sum_{i,t} \# \{ | e_{it} | > t \} } which is computationally conducted
 by repeatedly applying the median polish procedure in which one cell is
omitted (Davies, 2012; Terbeck & Davies, 1998).

Effect sizes of invariance are calculated as R-squared measures
of explained item slopes and intercepts after linking
in comparison to item parameters across groups
(Asparouhov & Muthen, 2014).

The function \eqn{linking.haberman.lq} uses the loss function \eqn{\rho(x)=|x|^q}.
The originally proposed Haberman linking can be obtained with \code{pow=2} (\eqn{q=2}).
The powers can also be estimated (argument \code{est_pow=TRUE}).
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\value{
A list with following entries
  \item{transf.pars}{Data frame with transformation parameters
\eqn{A_t} and \eqn{B_t}
  }
  \item{transf.personpars}{Data frame with linear transformation functions
  for person parameters}
  \item{joint.itempars}{Estimated joint item parameters \eqn{a_i} and \eqn{b_i}
      }
  \item{a.trans}{Transformed \eqn{a_{it}} parameters}
  \item{b.trans}{Transformed \eqn{b_{it}} parameters}
  \item{a.orig}{Original \eqn{a_{it}} parameters}
  \item{b.orig}{Original \eqn{b_{it}} parameters}
  \item{a.resid}{Residual \eqn{a_{it}} parameters (DIF parameters)}
  \item{b.resid}{Residual \eqn{b_{it}} parameters (DIF parameters)}
  \item{personpars}{Transformed person parameters}
  \item{es.invariance}{Effect size measures of invariance,
        separately for item slopes and intercepts.
        In the rows, \eqn{R^2} and \eqn{\sqrt{1-R^2}} are reported.}
  \item{es.robust}{Effect size measures of invariance based on
  robust estimation (if used).}
  \item{selitems}{Indices of items which are present in more than one
  study.}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\references{
Asparouhov, T., & Muthen, B. (2014). Multiple-group factor analysis alignment.
\emph{Structural Equation Modeling, 21}(4), 1-14.
\doi{10.1080/10705511.2014.919210}

Battauz, M. (2017). Multiple equating of separate IRT calibrations.
\emph{Psychometrika, 82}(3), 610-636.
\doi{10.1007/s11336-016-9517-x}

Davies, P. L. (2012). Interactions in the analysis of variance.
\emph{Journal of the American Statistical Association, 107}(500), 1502-1509.
\doi{10.1080/01621459.2012.726895}

Fox, J. (2015). \emph{Applied regression analysis and generalized linear models}.
Thousand Oaks: Sage.

Haberman, S. J. (2009). \emph{Linking parameter estimates derived
from an item response model through separate calibrations}.
ETS Research Report ETS RR-09-40. Princeton, ETS.
\doi{10.1002/j.2333-8504.2009.tb02197.x}

Kolen, M. J., & Brennan, R. L. (2014). \emph{Test equating, scaling, and linking:
Methods and practices}. New York: Springer.
\doi{10.1007/978-1-4939-0317-7}

Magis, D., & De Boeck, P. (2012). A robust outlier approach to prevent type I error
inflation in differential item functioning.
\emph{Educational and Psychological Measurement, 72}(2), 291-311.
\doi{10.1177/0013164411416975}

Maronna, R. A., Martin, R. D., & Yohai, V. J. (2006). \emph{Robust statistics}.
West Sussex: Wiley. \doi{10.1002/0470010940}

Terbeck, W., & Davies, P. L. (1998). Interactions and outliers in the two-way
analysis of variance. \emph{Annals of Statistics, 26}(4), 1279-1305.
doi: 10.1214/aos/1024691243

Tukey, J. W. (1977). \emph{Exploratory data analysis}. Addison-Wesley.

Weeks, J. P. (2010). \pkg{plink}: An \R package for linking mixed-format tests
using IRT-based methods. \code{Journal of Statistical Software, 35}(12), 1-33.
\doi{10.18637/jss.v035.i12}
}


%\author{
%Alexander Robitzsch
%}



%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See the \pkg{plink} package (Weeks, 2010) for a diversity of linking methods.

Mean-mean linking, Stocking-Lord and Haebara linking (see Kolen & Brennan, 2014,
for an overview) in the generalized logistic item response model can be conducted with
\code{\link{equating.rasch}}. See also \code{\link[TAM:tam.linking]{TAM::tam.linking}}
in the \pkg{TAM} package. Haebara linking and a robustified version of it can be
found in \code{\link{linking.haebara}}.

The invariance alignment method employs an optimization function based on
pairwise loss functions of item parameters (Asparouhov & Muthen, 2014),
see \code{\link{invariance.alignment}}.
}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\examples{
#############################################################################
# EXAMPLE 1: Item parameters data.pars1.rasch and data.pars1.2pl
#############################################################################

# Model 1: Linking three studies calibrated by the Rasch model
data(data.pars1.rasch)
mod1 <- sirt::linking.haberman( itempars=data.pars1.rasch )
summary(mod1)

# Model 1b: Linking these studies but weigh these studies by
#     proportion weights 3 : 0.5 : 1 (see below).
#     All weights are the same for each item but they could also
#     be item specific.
itempars <- data.pars1.rasch
itempars$wgt <- 1
itempars[ itempars$study=="study1","wgt"] <- 3
itempars[ itempars$study=="study2","wgt"] <- .5
mod1b <- sirt::linking.haberman( itempars=itempars )
summary(mod1b)

# Model 2: Linking three studies calibrated by the 2PL model
data(data.pars1.2pl)
mod2 <- sirt::linking.haberman( itempars=data.pars1.2pl )
summary(mod2)

# additive model instead of logarithmic model for item slopes
mod2b <- sirt::linking.haberman( itempars=data.pars1.2pl, a_log=FALSE )
summary(mod2b)

\dontrun{
#############################################################################
# EXAMPLE 2: Linking longitudinal data
#############################################################################
data(data.long)

#******
# Model 1: Scaling with the 1PL model

# scaling at T1
dat1 <- data.long[, grep("T1", colnames(data.long) ) ]
resT1 <- sirt::rasch.mml2( dat1 )
itempartable1 <- data.frame( "study"="T1", resT1$item[, c("item", "a", "b" ) ] )
# scaling at T2
dat2 <- data.long[, grep("T2", colnames(data.long) ) ]
resT2 <- sirt::rasch.mml2( dat2 )
summary(resT2)
itempartable2 <- data.frame( "study"="T2", resT2$item[, c("item", "a", "b" ) ] )
itempartable <- rbind( itempartable1, itempartable2 )
itempartable[,2] <- substring( itempartable[,2], 1, 2 )
# estimate linking parameters
mod1 <- sirt::linking.haberman( itempars=itempartable )

#******
# Model 2: Scaling with the 2PL model

# scaling at T1
dat1 <- data.long[, grep("T1", colnames(data.long) ) ]
resT1 <- sirt::rasch.mml2( dat1, est.a=1:6)
itempartable1 <- data.frame( "study"="T1", resT1$item[, c("item", "a", "b" ) ] )

# scaling at T2
dat2 <- data.long[, grep("T2", colnames(data.long) ) ]
resT2 <- sirt::rasch.mml2( dat2, est.a=1:6)
summary(resT2)
itempartable2 <- data.frame( "study"="T2", resT2$item[, c("item", "a", "b" ) ] )
itempartable <- rbind( itempartable1, itempartable2 )
itempartable[,2] <- substring( itempartable[,2], 1, 2 )
# estimate linking parameters
mod2 <- sirt::linking.haberman( itempars=itempartable )

#############################################################################
# EXAMPLE 3: 2 Studies - 1PL and 2PL linking
#############################################################################
set.seed(789)
I <- 20        # number of items
N <- 2000       # number of persons
# define item parameters
b <- seq( -1.5, 1.5, length=I )
# simulate data
dat1 <- sirt::sim.raschtype( stats::rnorm( N, mean=0,sd=1 ), b=b )
dat2 <- sirt::sim.raschtype( stats::rnorm( N, mean=0.5,sd=1.50 ), b=b )

#*** Model 1: 1PL
# 1PL Study 1
mod1 <- sirt::rasch.mml2( dat1, est.a=rep(1,I) )
summary(mod1)
# 1PL Study 2
mod2 <- sirt::rasch.mml2( dat2, est.a=rep(1,I) )
summary(mod2)

# collect item parameters
dfr1 <- data.frame( "study1", mod1$item$item, mod1$item$a, mod1$item$b )
dfr2 <- data.frame( "study2", mod2$item$item, mod2$item$a, mod2$item$b )
colnames(dfr2) <- colnames(dfr1) <- c("study", "item", "a", "b" )
itempars <- rbind( dfr1, dfr2 )

# Haberman linking
linkhab1 <- sirt::linking.haberman(itempars=itempars)
  ## Transformation parameters (Haberman linking)
  ##    study    At     Bt
  ## 1 study1 1.000  0.000
  ## 2 study2 1.465 -0.512
  ##
  ## Linear transformation for item parameters a and b
  ##    study   A_a   A_b    B_b
  ## 1 study1 1.000 1.000  0.000
  ## 2 study2 0.682 1.465 -0.512
  ##
  ## Linear transformation for person parameters theta
  ##    study A_theta B_theta
  ## 1 study1   1.000   0.000
  ## 2 study2   1.465   0.512
  ##
  ## R-Squared Measures of Invariance
  ##        slopes intercepts
  ## R2          1     0.9979
  ## sqrtU2      0     0.0456

#*** Model 2: 2PL
# 2PL Study 1
mod1 <- sirt::rasch.mml2( dat1, est.a=1:I )
summary(mod1)
# 2PL Study 2
mod2 <- sirt::rasch.mml2( dat2, est.a=1:I )
summary(mod2)

# collect item parameters
dfr1 <- data.frame( "study1", mod1$item$item, mod1$item$a, mod1$item$b )
dfr2 <- data.frame( "study2", mod2$item$item, mod2$item$a, mod2$item$b )
colnames(dfr2) <- colnames(dfr1) <- c("study", "item", "a", "b" )
itempars <- rbind( dfr1, dfr2 )

# Haberman linking
linkhab2 <- sirt::linking.haberman(itempars=itempars)
  ## Transformation parameters (Haberman linking)
  ##    study    At     Bt
  ## 1 study1 1.000  0.000
  ## 2 study2 1.468 -0.515
  ##
  ## Linear transformation for item parameters a and b
  ##    study   A_a   A_b    B_b
  ## 1 study1 1.000 1.000  0.000
  ## 2 study2 0.681 1.468 -0.515
  ##
  ## Linear transformation for person parameters theta
  ##    study A_theta B_theta
  ## 1 study1   1.000   0.000
  ## 2 study2   1.468   0.515
  ##
  ## R-Squared Measures of Invariance
  ##        slopes intercepts
  ## R2     0.9984     0.9980
  ## sqrtU2 0.0397     0.0443

#############################################################################
# EXAMPLE 4: 3 Studies - 1PL and 2PL linking
#############################################################################
set.seed(789)
I <- 20         # number of items
N <- 1500       # number of persons
# define item parameters
b <- seq( -1.5, 1.5, length=I )
# simulate data
dat1 <- sirt::sim.raschtype( stats::rnorm( N, mean=0, sd=1), b=b )
dat2 <- sirt::sim.raschtype( stats::rnorm( N, mean=0.5, sd=1.50), b=b )
dat3 <- sirt::sim.raschtype( stats::rnorm( N, mean=-0.2, sd=0.8), b=b )
# set some items to non-administered
dat3 <- dat3[, -c(1,4) ]
dat2 <- dat2[, -c(1,2,3) ]

#*** Model 1: 1PL in sirt
# 1PL Study 1
mod1 <- sirt::rasch.mml2( dat1, est.a=rep(1,ncol(dat1)) )
summary(mod1)
# 1PL Study 2
mod2 <- sirt::rasch.mml2( dat2, est.a=rep(1,ncol(dat2)) )
summary(mod2)
# 1PL Study 3
mod3 <- sirt::rasch.mml2( dat3, est.a=rep(1,ncol(dat3)) )
summary(mod3)

# collect item parameters
dfr1 <- data.frame( "study1", mod1$item$item, mod1$item$a, mod1$item$b )
dfr2 <- data.frame( "study2", mod2$item$item, mod2$item$a, mod2$item$b )
dfr3 <- data.frame( "study3", mod3$item$item, mod3$item$a, mod3$item$b )
colnames(dfr3) <- colnames(dfr2) <- colnames(dfr1) <- c("study", "item", "a", "b" )
itempars <- rbind( dfr1, dfr2, dfr3 )

# use person parameters
personpars <- list( mod1$person[, c("EAP","SE.EAP") ], mod2$person[, c("EAP","SE.EAP") ],
    mod3$person[, c("EAP","SE.EAP") ] )

# Haberman linking
linkhab1 <- sirt::linking.haberman(itempars=itempars, personpars=personpars)
# compare item parameters
round( cbind( linkhab1$joint.itempars[,-1], linkhab1$b.trans )[1:5,], 3 )
  ##            aj     bj study1 study2 study3
  ##   I0001 0.998 -1.427 -1.427     NA     NA
  ##   I0002 0.998 -1.290 -1.324     NA -1.256
  ##   I0003 0.998 -1.140 -1.068     NA -1.212
  ##   I0004 0.998 -0.986 -1.003 -0.969     NA
  ##   I0005 0.998 -0.869 -0.809 -0.872 -0.926

# summary of person parameters of second study
round( psych::describe( linkhab1$personpars[[2]] ), 2 )
  ##   var    n mean   sd median trimmed  mad   min  max range  skew kurtosis
  ## EAP      1 1500 0.45 1.36   0.41    0.47 1.52 -2.61 3.25  5.86 -0.08    -0.62
  ## SE.EAP   2 1500 0.57 0.09   0.53    0.56 0.04  0.49 0.84  0.35  1.47     1.56
  ##          se
  ## EAP    0.04
  ## SE.EAP 0.00

#*** Model 2: 2PL in TAM
library(TAM)
# 2PL Study 1
mod1 <- TAM::tam.mml.2pl( resp=dat1, irtmodel="2PL" )
pvmod1 <- TAM::tam.pv(mod1, ntheta=300, normal.approx=TRUE) # draw plausible values
summary(mod1)
# 2PL Study 2
mod2 <- TAM::tam.mml.2pl( resp=dat2, irtmodel="2PL" )
pvmod2 <- TAM::tam.pv(mod2, ntheta=300, normal.approx=TRUE)
summary(mod2)
# 2PL Study 3
mod3 <- TAM::tam.mml.2pl( resp=dat3, irtmodel="2PL" )
pvmod3 <- TAM::tam.pv(mod3, ntheta=300, normal.approx=TRUE)
summary(mod3)

# collect item parameters
#!!  Note that in TAM the parametrization is a*theta - b while linking.haberman
#!!  needs the parametrization a*(theta-b)
dfr1 <- data.frame( "study1", mod1$item$item, mod1$B[,2,1], mod1$xsi$xsi / mod1$B[,2,1] )
dfr2 <- data.frame( "study2", mod2$item$item, mod2$B[,2,1], mod2$xsi$xsi / mod2$B[,2,1] )
dfr3 <- data.frame( "study3", mod3$item$item, mod3$B[,2,1], mod3$xsi$xsi / mod3$B[,2,1] )
colnames(dfr3) <- colnames(dfr2) <- colnames(dfr1) <- c("study", "item", "a", "b" )
itempars <- rbind( dfr1, dfr2, dfr3 )

# define list containing person parameters
personpars <- list(  pvmod1$pv[,-1], pvmod2$pv[,-1], pvmod3$pv[,-1] )

# Haberman linking
linkhab2 <- sirt::linking.haberman(itempars=itempars,personpars=personpars)
  ##   Linear transformation for person parameters theta
  ##      study A_theta B_theta
  ##   1 study1   1.000   0.000
  ##   2 study2   1.485   0.465
  ##   3 study3   0.786  -0.192

# extract transformed person parameters
personpars.trans <- linkhab2$personpars

#############################################################################
# EXAMPLE 5: Linking with simulated item parameters containing outliers
#############################################################################

# simulate some parameters
I <- 38
set.seed(18785)
b <- stats::rnorm( I, mean=.3, sd=1.4 )
# simulate DIF effects plus some outliers
bdif <- stats::rnorm(I,mean=.4,sd=.09)+( stats::runif(I)>.9 )* rep( 1*c(-1,1)+.4, each=I/2 )
# create item parameter table
itempars <- data.frame( "study"=paste0("study",rep(1:2, each=I)),
                "item"=paste0( "I", 100 + rep(1:I,2) ), "a"=1,
                 "b"=c( b, b + bdif  )  )

#*** Model 1: Haberman linking with least squares regression
mod1 <- sirt::linking.haberman( itempars=itempars )
summary(mod1)

#*** Model 2: Haberman linking with robust bisquare regression with fixed trimming value
mod2 <- sirt::linking.haberman( itempars=itempars, estimation="BSQ", b_trim=.4)
summary(mod2)

#*** Model 2: Haberman linking with robust bisquare regression with estimated trimming value
mod3 <- sirt::linking.haberman( itempars=itempars, estimation="BSQ")
summary(mod3)

## see also Example 3 of ?sirt::robust.linking

#############################################################################
# EXAMPLE 6: Toy example of Magis and De Boeck (2012)
#############################################################################

# define item parameters from Magis & De Boeck (20212, p. 293)
b1 <- c(1,1,1,1)
b2 <- c(1,1,1,2)
itempars <- data.frame(study=rep(1:2, each=4), item=rep(1:4,2), a=1, b=c(b1,b2) )

#- Least squares regression
mod1 <- sirt::linking.haberman( itempars=itempars, estimation="OLS")
summary(mod1)

#- Bisquare regression with estimated and fixed trimming factors
mod2 <- sirt::linking.haberman( itempars=itempars, estimation="BSQ")
mod2a <- sirt::linking.haberman( itempars=itempars, estimation="BSQ", b_trim=.4)
mod2b <- sirt::linking.haberman( itempars=itempars, estimation="BSQ", b_trim=1.2)
summary(mod2)
summary(mod2a)
summary(mod2b)

#- Least squares trimmed regression
mod3 <- sirt::linking.haberman( itempars=itempars, estimation="LTS")
summary(mod3)

#- median regression
mod4 <- sirt::linking.haberman( itempars=itempars, estimation="MED")
summary(mod4)

#############################################################################
# EXAMPLE 7: Simulated example with directional DIF
#############################################################################

set.seed(98)
I <- 8
mu <- c(-.5, 0, .5)
b <- sample(seq(-1.5,1.5, len=I))
sd_dif <- 0.001
pars <- outer(b, mu, "+") + stats::rnorm(I*3, sd=sd_dif)
ind <- c(1,2); pars[ind,1] <- pars[ind,1] + c(.5,.5)
ind <- c(3,4); pars[ind,2] <- pars[ind,2] + (-1)*c(.6,.6)
ind <- c(5,6); pars[ind,3] <- pars[ind,3] + (-1)*c(1,1)

# median polish (=stats::medpolish())
tmod1 <- sirt:::L1_polish(x=pars)
# L0 polish with tolerance criterion of .3
tmod2 <- sirt::L0_polish(x=pars, tol=.3)

#- prepare itempars input
itempars <- sirt::linking_haberman_itempars_prepare(b=pars)

#- compare different estimation functions for Haberman linking
mod01 <- sirt::linking.haberman(itempars, estimation="L1")
mod02 <- sirt::linking.haberman(itempars, estimation="L0", b_trim=.3)
mod1 <- sirt::linking.haberman(itempars, estimation="OLS")
mod2 <- sirt::linking.haberman(itempars, estimation="BSQ")
mod2a <- sirt::linking.haberman(itempars, estimation="BSQ", b_trim=.4)
mod3 <- sirt::linking.haberman(itempars, estimation="MED")
mod4 <- sirt::linking.haberman(itempars, estimation="LTS")
mod5 <- sirt::linking.haberman(itempars, estimation="HUB")
mod01$transf.pars
mod02$transf.pars
mod1$transf.pars
mod2$transf.pars
mod2a$transf.pars
mod3$transf.pars
mod4$transf.pars
mod5$transf.pars

#############################################################################
# EXAMPLE 8: Many studies and directional DIF
#############################################################################

## dataset 2
set.seed(98)
I <- 10 # number of items
S <- 7  # number of studies
mu <- round( seq(0, 1, len=S))
b <- sample(seq(-1.5,1.5, len=I))
sd_dif <- 0.001
pars0 <- pars <- outer(b, mu, "+") + stats::rnorm(I*S, sd=sd_dif)

# select n_dif items at random per group and set it to dif or -dif
n_dif <- 2
dif <- .6
for (ss in 1:S){
    ind <- sample( 1:I, n_dif )
    pars[ind,ss] <- pars[ind,ss] + dif*sign( runif(1) - .5 )
}

# check DIF
pars - pars0

#* estimate models
itempars <- sirt::linking_haberman_itempars_prepare(b=pars)
mod0 <- sirt::linking.haberman(itempars, estimation="L0", b_trim=.2)
mod1 <- sirt::linking.haberman(itempars, estimation="OLS")
mod2 <- sirt::linking.haberman(itempars, estimation="BSQ")
mod2a <- sirt::linking.haberman(itempars, estimation="BSQ", b_trim=.4)
mod3 <- sirt::linking.haberman(itempars, estimation="MED")
mod3a <- sirt::linking.haberman(itempars, estimation="L1")
mod4 <- sirt::linking.haberman(itempars, estimation="LTS")
mod5 <- sirt::linking.haberman(itempars, estimation="HUB")
mod0$transf.pars
mod1$transf.pars
mod2$transf.pars
mod2a$transf.pars
mod3$transf.pars
mod3a$transf.pars
mod4$transf.pars
mod5$transf.pars

#* compare results with Haebara linking
mod11 <- sirt::linking.haebara(itempars, dist="L2")
mod12 <- sirt::linking.haebara(itempars, dist="L1")
summary(mod11)
summary(mod12)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Linking}
%% \keyword{Equating}
