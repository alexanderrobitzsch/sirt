%% File Name: invariance.alignment.Rd
%% File Version: 1.39

\name{invariance.alignment}
\alias{invariance.alignment}
\alias{summary.invariance.alignment}
\alias{plot.invariance.alignment}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Alignment Procedure for Linking under Approximate Invariance
}
\description{
This function does alignment under approximate
invariance for \eqn{G} groups and \eqn{I} items
(Asparouhov & Muthen, 2014; Muthen, 2014).
It is assumed that item loadings and intercepts are
previously estimated under the assumption of a factor
with zero mean and a variance of one. 
}

%% \emph{Note that this function is in an early experimental stage.}


\usage{
invariance.alignment(lambda, nu, wgt=NULL, align.scale = c(1, 1), 
    align.pow = c(1, 1), eps = .01, h = 0.001, max.increment = 0.2, 
    increment.factor = c(1.001,1.02,1.04,1.08), maxiter = 300, conv = 1e-04, 
    fac.oldpar= c(.01,.2,.5,.85), psi0.init = NULL, alpha0.init = NULL, progress = TRUE)

\method{summary}{invariance.alignment}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda}{
A \eqn{G \times I} matrix with item loadings
}
  \item{nu}{
A \eqn{G \times I}  matrix with item intercepts
}
  \item{wgt}{
A \eqn{G \times I} matrix for weighing groups 
for each item
}
  \item{align.scale}{
A vector of length two containing scale parameter
\eqn{a_\lambda} and \eqn{a_\nu} (see Details) 
}
  \item{align.pow}{
A vector of length two containing power
\eqn{p_\lambda} and \eqn{p_\nu} (see Details)
}
  \item{eps}{
A parameter in the optimization function
}
  \item{h}{
Numerical differentiation parameter 
}
  \item{max.increment}{
Maximum increment in each iteration
}
  \item{increment.factor}{
A numerical larger than one indicating the extent
of the decrease of \code{max.increment} in every
iteration.
}
  \item{maxiter}{
Maximum number of iterations
}
  \item{conv}{
Maximum parameter change of the optimization function
}
\item{fac.oldpar}{Convergence acceleration parameter between 0 and 1. 
    This parameter defines the relative weight the previous parameter value
    for the calculation of the parameter update.
    The default is .85. But experiment with this value and study
    the obtained results.
}
  \item{psi0.init}{
An optional vector of initial \eqn{\psi_0} parameters
}
  \item{alpha0.init}{
An optional vector of initial \eqn{\alpha_0} parameters
}
  \item{progress}{
An optional logical indicating whether computational
should be printed.
}
  \item{object}{
Object of class \code{invariance.alignment}
  }
  \item{\dots}{
Further optional arguments to be passed}
}

\details{
For \eqn{G} groups and \eqn{I} items, item loadings \eqn{\lambda_{ig0}}
and intercepts \eqn{\nu_{ig0}} are available and have been estimated
in a 1-dimensional factor analysis assuming a standardized factor.

The alignment procedure searches means \eqn{\alpha_{g0}}
and standard deviations \eqn{\psi_{g0}} using an alignment
optimization function \eqn{F}. This function is defined as
\deqn{F = \sum_i \sum_{ g_1 < g_2} w_{i,g1} w_{i,g2} 
    f_\lambda( \lambda_{i g_1 ,1} - \lambda_{i g_2 ,1} )
    + \sum_i \sum_{ g_1 < g_2} w_{i,g1} w_{i,g2} 
f_\nu( \nu_{i g_1 ,1} - \nu_{i g_2 ,1} ) }
where the aligned item parameters \eqn{\lambda_{i g ,1}}
and \eqn{\nu_{i g ,1}} are defined such that

\deqn{ \lambda_{i g ,1} = \lambda_{i g 0} / \psi_{g0} 
    \qquad \mbox{and} \qquad 
    \nu_{i g ,1} = \nu_{i g 0} -  \alpha_{g0} \lambda_{ig0} / \psi_{g0}
            }
and the optimization functions are defined as
\deqn{ f_\lambda (x) = [ ( x/ a_\lambda )^2 + \varepsilon ]^{p_\lambda}
    \qquad \mbox{and} \qquad
    f_\nu (x) = [ ( x/ a_\nu )^2 + \varepsilon ]^{p_\nu}
            }
using a small \eqn{ \varepsilon > 0} (e.g. .0001) to obtain
a differentiable optimization function. 

For identification reasons, the product \eqn{\Pi_g \psi_{g0}}
of all group standard deviations is set to one. The mean
\eqn{\alpha_{g0}} of the first group is set to zero.

Note that the standard deviations \eqn{\psi_g} are estimated
due to minimizing the sum of \eqn{f_\lambda} functions while
means \eqn{\alpha_g} are obtained by minimizing the \eqn{f_\nu}
part with fixed \eqn{\psi_g} parameters. Therefore, the original approach
of Asparouhov and Muthen (2014) is split into a two-step procedure.

Note that Asparouhov and Muthen (2014) use \eqn{a_\lambda=a_\nu=1}
(which can be modified in \code{align.scale})
and \eqn{p_\lambda = p_\nu = 1/4} (which can be modified in \code{align.pow}).
In case of \eqn{p_\lambda=1}, the penalty is approximately 
\eqn{f_\lambda(x) = x^2 }, in case of \eqn{p_\lambda=1/4}
it is approximately \eqn{f_\lambda(x) = \sqrt{|x|} }.

Effect sizes of approximate invariance based on \eqn{R^2} have
been proposed by Asparouhov and Muthen (2014). These are
calculated separately for item loading and intercepts, resulting
in \eqn{R^2_\lambda} and \eqn{R^2_\nu} measures which are
included in the output \code{es.invariance}. In addition, 
the average correlation of aligned item parameters among groups (\code{rbar})
is reported.

\emph{Metric invariance} means that all aligned item loadings \eqn{\lambda_{ig,1}}
are equal across groups and therefore \eqn{R^2_\lambda = 1}.
\emph{Scalar invariance} means that all aligned item loadings 
\eqn{\lambda_{ig,1}} and aligned item intercepts \eqn{\nu_{ig,1}} are
equal across groups and therefore \eqn{R^2_\lambda = 1} and \eqn{R^2_\nu = 1}
(see Vandenberg & Lance, 2000).
}

\value{
A list with following entries
\item{pars}{Aligned distribution parameters}
\item{itempars.aligned}{Aligned item parameters for all groups}
\item{es.invariance}{Effect sizes of approximate invariance}
\item{lambda.aligned}{Aligned \eqn{ \lambda_{i g ,1}} parameters}
\item{lambda.resid}{Residuals of \eqn{ \lambda_{i g ,1}} parameters}
\item{nu.aligned}{Aligned \eqn{ \nu_{i g ,1}} parameters}
\item{nu.resid}{Residuals of \eqn{ \nu_{i g ,1}} parameters}
\item{Niter}{Number of iterations for \eqn{f_\lambda} and
\eqn{f_\nu} optimization functions}
\item{miniter}{Iteration index with minimum optimization value}
\item{fopt}{Minimum optimization value}
\item{align.scale}{Used alignment scale parameters}
\item{align.pow}{Used alignment power parameters}
}


\references{
Asparouhov, T., & Muthen, B. (2014).
Multiple-group factor analysis alignment.
\emph{Structural Equation Modeling}, \bold{21}, 1-14.
\url{http://www.statmodel.com/Alignment.shtml}

Muthen, B., & Asparouhov, T. (2014). IRT studies of many groups: The alignment method. 
\emph{Frontiers in Psychology | Quantitative Psychology and Measurement},
\bold{5:978}. doi: 10.3389/fpsyg.2014.00978,
\url{http://journal.frontiersin.org/Journal/10.3389/fpsyg.2014.00978/abstract}.

Vandenberg, R. J., & Lance, C. E. (2000). A review and synthesis of the
measurement invariance literature: Suggestions, practices, and
recommendations for organizational research. \emph{Organizational Research
Methods}, \bold{3}, 4-70.
}

\author{
Alexander Robitzsch
}

%\note{
%xxxx
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For IRT linking see also \code{\link{linking.haberman}}.

For modeling random item effects for loadings and intercepts
see \code{\link{mcmc.2pno.ml}}.
}

\examples{
#############################################################################
# EXAMPLE 1: Item parameters cultural activities 
#############################################################################

data( data.activity.itempars )
lambda <- data.activity.itempars$lambda
nu <- data.activity.itempars$nu
Ng <-  data.activity.itempars$N
wgt <- matrix( sqrt(Ng) , length(Ng) , ncol(nu) )

#***
# Model 1: Alignment using a quadratic loss function
#   -> use the default of align.pow=c(1,1) and align.scale=c(1,1)
mod1 <- sirt::invariance.alignment( lambda , nu , wgt )
summary(mod1)
  ##   Effect Sizes of Approximate Invariance 
  ##          loadings intercepts
  ##   R2       0.9944     0.9988
  ##   sqrtU2   0.0748     0.0346
  ##   rbar     0.9265     0.9735

#****
# Model 2: Different powers for alignment
mod2 <- sirt::invariance.alignment( lambda , nu , wgt ,  align.pow=c(.25,1/2) , 
              align.scale=c(.95,.95) ,  max.increment=.1)
summary(mod2)

# compare means from Models 1 and 2
plot( mod1$pars$alpha0 , mod2$pars$alpha0 , pch=16 ,
    xlab= "M (Model 1)" , ylab="M (Model 2)" , xlim=c(-.3,.3) , ylim=c(-.3,.3) )
lines( c(-1,1) , c(-1,1) , col="gray")    
round( cbind( mod1$pars$alpha0 , mod2$pars$alpha0 ) , 3 )
round( mod1$nu.resid , 3)
round( mod2$nu.resid ,3 )

#****
# Model 3: Low powers for alignment of scale and power
# Note that setting increment.factor larger than 1 seems necessary
mod3 <- sirt::invariance.alignment( lambda , nu , wgt , align.pow=c(.25,.35) , 
            align.scale=c(.55,.55) , psi0.init=mod1$psi0 , alpha0.init = mod1$alpha0 )
summary(mod3)

# compare mean and SD estimates of Models 1 and 3
plot( mod1$pars$alpha0 , mod3$pars$alpha0 , pch=16)
plot( mod1$pars$psi0 , mod3$pars$psi0 , pch=16)

# compare residuals for Models 1 and 3   
# plot lambda
plot( abs(as.vector(mod1$lambda.resid)) , abs(as.vector(mod3$lambda.resid)) ,
      pch=16 , xlab="Residuals lambda (Model 1)" ,
      ylab="Residuals lambda (Model 3)" , xlim=c(0,.1) , ylim=c(0,.1))
lines( c(-3,3),c(-3,3) , col="gray")    
# plot nu
plot( abs(as.vector(mod1$nu.resid)) , abs(as.vector(mod3$nu.resid)) ,
      pch=16 , xlab="Residuals nu (Model 1)" , ylab="Residuals nu (Model 3)" , 
      xlim=c(0,.4),ylim=c(0,.4))
lines( c(-3,3),c(-3,3) , col="gray")

\dontrun{
#############################################################################
# EXAMPLE 2: Comparison 4 groups | data.inv4gr 
#############################################################################

data(data.inv4gr)
dat <- data.inv4gr
miceadds::library_install("semTools")

model1 <- "
    F =~ I01 + I02 + I03 + I04 + I05 + I06 + I07 + I08 + I09 + I10 + I11
    F ~~ 1*F
    "

res <- semTools::measurementInvariance(model1, std.lv =TRUE , data=dat , group="group")
  ##   Measurement invariance tests:
  ##   
  ##   Model 1: configural invariance:
  ##       chisq        df    pvalue       cfi     rmsea       bic 
  ##     162.084   176.000     0.766     1.000     0.000 95428.025 
  ##   
  ##   Model 2: weak invariance (equal loadings):
  ##       chisq        df    pvalue       cfi     rmsea       bic 
  ##     519.598   209.000     0.000     0.973     0.039 95511.835 
  ##   
  ##   [Model 1 versus model 2]
  ##     delta.chisq      delta.df delta.p.value     delta.cfi 
  ##         357.514        33.000         0.000         0.027 
  ##   
  ##   Model 3: strong invariance (equal loadings + intercepts):
  ##       chisq        df    pvalue       cfi     rmsea       bic 
  ##    2197.260   239.000     0.000     0.828     0.091 96940.676 
  ##   
  ##   [Model 1 versus model 3]
  ##     delta.chisq      delta.df delta.p.value     delta.cfi 
  ##        2035.176        63.000         0.000         0.172 
  ##   
  ##   [Model 2 versus model 3]
  ##     delta.chisq      delta.df delta.p.value     delta.cfi 
  ##        1677.662        30.000         0.000         0.144 
  ##   

# extract item parameters separate group analyses
ipars <- lavaan::parameterEstimates(res$fit.configural)
# extract lambda's: groups are in rows, items in columns
lambda <- matrix( ipars[ ipars$op == "=~" , "est"] , nrow=4 ,  byrow=TRUE)
colnames(lambda) <- colnames(dat)[-1]
# extract nu's
nu <- matrix( ipars[ ipars$op == "~1"  & ipars$se != 0 , "est" ], nrow=4 ,  byrow=TRUE)
colnames(nu) <- colnames(dat)[-1]

# Model 1: least squares optimization
mod1 <- sirt::invariance.alignment( lambda=lambda , nu=nu )
summary(mod1)
  ##   Effect Sizes of Approximate Invariance 
  ##          loadings intercepts
  ##   R2       0.9826     0.9972
  ##   sqrtU2   0.1319     0.0526
  ##   rbar     0.6237     0.7821
  ##   -----------------------------------------------------------------
  ##   Group Means and Standard Deviations 
  ##     alpha0  psi0
  ##   1  0.000 0.965
  ##   2 -0.105 1.098
  ##   3 -0.081 1.011
  ##   4  0.171 0.935

# Model 2: sparse target function
mod2 <- sirt::invariance.alignment( lambda=lambda , nu=nu , align.pow=c(1/4,1/4) )
summary(mod2)
  ##   Effect Sizes of Approximate Invariance 
  ##          loadings intercepts
  ##   R2       0.9824     0.9972
  ##   sqrtU2   0.1327     0.0529
  ##   rbar     0.6237     0.7856
  ##   -----------------------------------------------------------------
  ##   Group Means and Standard Deviations 
  ##     alpha0  psi0
  ##   1 -0.002 0.965
  ##   2 -0.107 1.098
  ##   3 -0.083 1.011
  ##   4  0.170 0.935

#############################################################################
# EXAMPLE 3: European Social Survey data.ess2005 
#############################################################################

data(data.ess2005)
lambda <- data.ess2005$lambda
nu <- data.ess2005$nu

# Model 1: least squares optimization
mod1 <- sirt::invariance.alignment( lambda=lambda , nu=nu )
summary(mod1)

# Model 2: sparse target function and definition of scales
mod2 <- sirt::invariance.alignment( lambda=lambda , nu=nu , align.pow=c(1/4,1/4) , 
            align.scale= c( .2 , .3)  )
summary(mod2)

# compare results of Model 1 and Model 2
round( cbind( mod1$pars , mod2$pars ) , 2 )
  ##      alpha0 psi0 alpha0 psi0
  ##   1    0.06 0.87   0.05 0.91
  ##   2   -0.51 1.03  -0.37 0.99
  ##   3    0.18 0.97   0.25 1.04
  ##   4   -0.67 0.90  -0.53 0.90
  ##   5    0.09 0.98   0.10 0.99
  ##   6    0.23 1.03   0.28 1.00
  ##   7    0.27 0.97   0.14 1.10
  ##   8    0.18 0.90   0.07 0.89
  ##   [...]

# look at nu residuals to explain differences in means
round( mod1$nu.resid , 2)
  ##         ipfrule ipmodst ipbhprp imptrad
  ##    [1,]    0.15   -0.25   -0.01    0.01
  ##    [2,]   -0.18    0.23    0.10   -0.24
  ##    [3,]    0.22   -0.34    0.05   -0.02
  ##    [4,]    0.29   -0.04    0.12   -0.53
  ##    [5,]   -0.32    0.19    0.00    0.13
  ##    [6,]    0.05   -0.21    0.05    0.04
  ##    [7,]   -0.26    0.54   -0.15   -0.02
  ##    [8,]    0.07   -0.05   -0.10    0.12
round( mod2$nu.resid , 2)
  ##         ipfrule ipmodst ipbhprp imptrad
  ##    [1,]    0.16   -0.25    0.00    0.02
  ##    [2,]   -0.27    0.14    0.00   -0.30
  ##    [3,]    0.18   -0.37    0.00   -0.05
  ##    [4,]    0.19   -0.13    0.00   -0.60
  ##    [5,]   -0.33    0.19   -0.01    0.12
  ##    [6,]    0.00   -0.23    0.00    0.01
  ##    [7,]   -0.16    0.64   -0.01    0.04
  ##    [8,]    0.15    0.02   -0.02    0.19

round( rowMeans( mod1$nu.resid )[1:8] , 2 )
  ##   [1] -0.02 -0.02 -0.02 -0.04  0.00 -0.02  0.03  0.01
round( rowMeans( mod2$nu.resid )[1:8] , 2 )
  ##   [1] -0.02 -0.11 -0.06 -0.14 -0.01 -0.06  0.13  0.09

#############################################################################
# EXAMPLE 4: Linking with item parameters containing outliers 
#############################################################################
    
# see Help file in linking.robust
    
# simulate some item difficulties in the Rasch model
I <- 38
set.seed(18785)
itempars <- data.frame("item" = paste0("I",1:I) )
itempars$study1 <- stats::rnorm( I , mean = .3 , sd =1.4 )
# simulate DIF effects plus some outliers
bdif <- stats::rnorm(I,mean=.4,sd=.09)+( stats::runif(I)>.9 )* rep( 1*c(-1,1)+.4 , each=I/2 )
itempars$study2 <- itempars$study1 + bdif
# create input for function invariance.alignment
nu <- t( itempars[,2:3] )
colnames(nu) <- itempars$item
lambda <- 1+0*nu

# linking using least squares optimization
mod1 <- sirt::invariance.alignment( lambda=lambda , nu=nu )
summary(mod1)
  ##   Group Means and Standard Deviations 
  ##          alpha0 psi0
  ##   study1 -0.286    1
  ##   study2  0.286    1

# linking using powers of .5
mod2 <- sirt::invariance.alignment( lambda=lambda , nu=nu , align.pow=c(.5,.5) )
summary(mod2)
  ##   Group Means and Standard Deviations 
  ##          alpha0 psi0
  ##   study1 -0.213    1
  ##   study2  0.213    1

# linking using powers of .25
mod3 <- sirt::invariance.alignment( lambda=lambda , nu=nu , align.pow=c(.25,.25) )
summary(mod3)
  ##   Group Means and Standard Deviations 
  ##          alpha0 psi0
  ##   study1 -0.207    1
  ##   study2  0.207    1

#############################################################################
# EXAMPLE 5: Linking gender groups with data.math 
#############################################################################

data(data.math)
dat <- data.math$data
dat.male <- dat[ dat$female == 0 , substring( colnames(dat) ,1,1) == "M"  ]
dat.female <- dat[ dat$female == 1 , substring( colnames(dat) ,1,1) == "M"  ]

#*************************
# Model 1: Linking using the Rasch model
mod1m <- sirt::rasch.mml2( dat.male )
mod1f <- sirt::rasch.mml2( dat.female )

# create objects for invariance.alignment
nu <- rbind( mod1m$item$thresh , mod1f$item$thresh )
colnames(nu) <- mod1m$item$item
rownames(nu) <- c("male" , "female")
lambda <- 1+0*nu

# mean of item difficulties
round( rowMeans(nu) , 3 )
  ##     male female 
  ##   -0.081 -0.049 

# Linking using least squares optimization
res1a <- sirt::invariance.alignment( lambda , nu , align.scale = c( .3 , .5 ) )
summary(res1a)
  ##   Effect Sizes of Approximate Invariance 
  ##          loadings intercepts
  ##   R2            1     0.9801
  ##   sqrtU2        0     0.1412
  ##   rbar          1     0.9626
  ##   -----------------------------------------------------------------
  ##   Group Means and Standard Deviations 
  ##          alpha0 psi0
  ##   male   -0.016    1
  ##   female  0.016    1

# Linking using optimization with absolute values
res1b <- sirt::invariance.alignment( lambda , nu , align.scale = c( .3 , .5 ) ,
                align.pow=c( .5 , .5 ) )
summary(res1b)
  ##   Group Means and Standard Deviations 
  ##          alpha0 psi0
  ##   male   -0.045    1
  ##   female  0.045    1

#-- compare results with Haberman linking
I <- ncol(dat.male)
itempartable <- data.frame( "study" = rep( c("male" , "female") , each=I ) )
itempartable$item <- c( paste0(mod1m$item$item)  ,  paste0(mod1f$item$item) )
itempartable$a <- 1
itempartable$b <- c( mod1m$item$b , mod1f$item$b )
# estimate linking parameters
res1c <- sirt::linking.haberman( itempars= itempartable )
  ##   Transformation parameters (Haberman linking)
  ##      study At     Bt
  ##   1 female  1  0.000
  ##   2   male  1 -0.032
  ##   Linear transformation for person parameters theta
  ##      study A_theta B_theta
  ##   1 female       1   0.000
  ##   2   male       1   0.032
  ##   R-Squared Measures of Invariance
  ##          slopes intercepts
  ##   R2          1     0.9801
  ##   sqrtU2      0     0.1412

#-- results of equating.rasch
x <- itempartable[ 1:I , c("item" , "b") ]
y <- itempartable[ I + 1:I , c("item" , "b") ]
res1d <- sirt::equating.rasch( x , y )
round( res1d$B.est , 3 )
  ##     Mean.Mean Haebara Stocking.Lord
  ##   1     0.032   0.032         0.029

#*************************
# Model 2: Linking using the 2PL model
I <- ncol(dat.male)
mod2m <- sirt::rasch.mml2( dat.male , est.a=1:I)
mod2f <- sirt::rasch.mml2( dat.female , est.a=1:I)

# create objects for invariance.alignment
nu <- rbind( mod2m$item$thresh , mod2f$item$thresh )
colnames(nu) <- mod2m$item$item
rownames(nu) <- c("male" , "female")
lambda <- rbind( mod2m$item$a , mod2f$item$a )
colnames(lambda) <- mod2m$item$item
rownames(lambda) <- c("male" , "female")

res2a <- sirt::invariance.alignment( lambda , nu , align.scale = c( .3 , .5 ) )
summary(res2a)
  ##   Effect Sizes of Approximate Invariance 
  ##          loadings intercepts
  ##   R2       0.9589     0.9682
  ##   sqrtU2   0.2027     0.1782
  ##   rbar     0.5177     0.9394
  ##   -----------------------------------------------------------------
  ##   Group Means and Standard Deviations 
  ##          alpha0  psi0
  ##   male   -0.044 0.968
  ##   female  0.047 1.034

res2b <- sirt::invariance.alignment( lambda , nu , align.scale = c( .3 , .5 ) ,
                align.pow=c( .5 , .5 ) )
summary(res2b)
  ##   Group Means and Standard Deviations 
  ##          alpha0  psi0
  ##   male   -0.046 1.053
  ##   female  0.041 0.951

# compare results with Haberman linking
I <- ncol(dat.male)
itempartable <- data.frame( "study" = rep( c("male" , "female") , each=I ) )
itempartable$item <- c( paste0(mod2m$item$item)  ,  paste0(mod2f$item$item ) )
itempartable$a <- c( mod2m$item$a , mod2f$item$a )
itempartable$b <- c( mod2m$item$b , mod2f$item$b )
# estimate linking parameters
res2c <- sirt::linking.haberman( itempars= itempartable )
  ##   Transformation parameters (Haberman linking)
  ##      study    At   Bt
  ##   1 female 1.000 0.00
  ##   2   male 1.041 0.09
  ##   Linear transformation for person parameters theta
  ##      study A_theta B_theta
  ##   1 female   1.000    0.00
  ##   2   male   1.041   -0.09
  ##   R-Squared Measures of Invariance
  ##          slopes intercepts
  ##   R2     0.9554     0.9484
  ##   sqrtU2 0.2111     0.2273
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Linking}
\keyword{Approximate Invariance}
\keyword{Alignment}
\keyword{summary}

