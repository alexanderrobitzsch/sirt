%% File Name: invariance.alignment.Rd
%% File Version: 1.712

\name{invariance.alignment}
\alias{invariance.alignment}
\alias{summary.invariance.alignment}
%% \alias{plot.invariance.alignment}
\alias{invariance_alignment_constraints}
\alias{summary.invariance_alignment_constraints}
\alias{invariance_alignment_simulate}
\alias{invariance_alignment_cfa_config}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Alignment Procedure for Linking under Approximate Invariance
}
\description{
The function \code{invariance.alignment} performs alignment under approximate
invariance for \eqn{G} groups and \eqn{I} items
(Asparouhov & Muthen, 2014; Byrne & van de Vijver, 2017; DeMars, 2020; Finch, 2016;
Fischer & Karl, 2019; Flake & McCoach, 2018; Kim et al., 2017; Marsh et al., 2018;
Muthen & Asparouhov, 2014, 2018; Pokropek, Davidov & Schmidt, 2019).
It is assumed that item loadings and intercepts are
previously estimated as a unidimensional factor model under the assumption of a factor
with zero mean and a variance of one.

The function \code{invariance_alignment_constraints} postprocesses the output of the
\code{invariance.alignment} function and estimates item parameters under equality
constraints for prespecified absolute values of parameter tolerance.

The function \code{invariance_alignment_simulate} simulates a one-factor model
for multiple groups for given matrices of \eqn{\nu} and \eqn{\lambda} parameters of
item intercepts and item slopes (see Example 6).

The function \code{invariance_alignment_cfa_config} estimates one-factor
models separately for each group as a preliminary step for invariance
alignment (see Example 6). Sampling weights are accommodated by the
argument \code{weights}.  The computed variance matrix \code{vcov} by this function
can be used to obtain standard errors in the \code{invariance.alignment} function
if it is supplied as the argument \code{vcov}.
}

%% \emph{Note that this function is in an early experimental stage.}


\usage{
invariance.alignment(lambda, nu, wgt=NULL, align.scale=c(1, 1),
    align.pow=c(.5, .5), eps=1e-3, psi0.init=NULL, alpha0.init=NULL, center=FALSE,
    optimizer="optim", fixed=NULL, meth=1, vcov=NULL, eps_grid=seq(0,-10, by=-.5),
    num_deriv=FALSE, ...)

\method{summary}{invariance.alignment}(object, digits=3, file=NULL, ...)

invariance_alignment_constraints(model, lambda_parm_tol, nu_parm_tol )

\method{summary}{invariance_alignment_constraints}(object, digits=3, file=NULL, ...)

invariance_alignment_simulate(nu, lambda, err_var, mu, sigma, N, output="data",
     groupwise=FALSE, exact=FALSE)

invariance_alignment_cfa_config(dat, group, weights=NULL, model="2PM", verbose=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda}{
A \eqn{G \times I} matrix with item loadings
}
  \item{nu}{
A \eqn{G \times I} matrix with item intercepts
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
  \item{psi0.init}{
An optional vector of initial \eqn{\psi_0} parameters
}
  \item{alpha0.init}{
An optional vector of initial \eqn{\alpha_0} parameters
}
 \item{center}{
Logical indicating whether estimated means and standard deviations should
be centered.
}
\item{optimizer}{Name of the optimizer chosen for alignment. Options are
   \code{"optim"} (using \code{\link[stats:optim]{stats::optim}})
   or \code{"nlminb"} (using \code{\link[stats:nlminb]{stats::nlminb}}).}
\item{fixed}{Logical indicating whether SD of first group should
be fixed to one. If \code{fixed=FALSE}, the product of all SDs is set to one.
If \code{NULL}, then \code{fixed} is automatically chosen by default. For many groups,
\code{fixed=FALSE} is chosen.}
\item{meth}{Type of method used for optimization function. \code{meth=1} is the default
and the optimization function used in Mplus. \code{meth=2} uses logarithmized
item loadings in alignment. The choice \code{meth=4} uses the constraint
\eqn{\prod_g \psi_g=1} and adds the penalty \eqn{\lambda \sum_g \alpha_g^2} for
a fixed value \eqn{\lambda} that depends on the weights \code{wgt}
(similar to Mplus' free method).
The choice \code{meth=3} only uses the constraint
\eqn{\prod_g \psi_g=1} (similar to Mplus' FIXED method).}
\item{vcov}{Variance matrix produced by \code{invariance_alignment_cfa_config}
for standard error computation. If a matrix is provided, standard errors
are computed.}
\item{eps_grid}{Grid of logarithmized epsilon values in optimization}
\item{num_deriv}{Logical indicating whether numerical derivatives should be used}
  \item{object}{
Object of class \code{invariance.alignment}
  }
\item{digits}{Number of digits used for rounding}
\item{file}{Optional file name in which summary should be sunk}
  \item{\dots}{
Further optional arguments to be passed}
\item{model}{Model of class \code{invariance.alignment}.
For \code{invariance_alignment_cfa_config}: Model type: \code{"2PM"} for two-parameter
model with unequal loadings
and \code{"1PM"} with equal loadings and equal residual variances}
\item{lambda_parm_tol}{Parameter tolerance for \eqn{\lambda} parameters}
\item{nu_parm_tol}{Parameter tolerance for \eqn{\nu} parameters}
\item{err_var}{Error variance}
\item{mu}{Vector of means}
\item{sigma}{Vector of standard deviations}
\item{N}{Vector of sample sizes per group}
\item{output}{Specifies output type: \code{"data"} for dataset and \code{"suffstat"}
for sufficient statistics (i.e., means and covariance matrices)}
\item{groupwise}{Logical indicating whether group-wise output is requested}
\item{exact}{Logical indicating whether distributions should be exactly preserved in
simulated data}
\item{dat}{Dataset with items or a list containing sufficient statistics}
\item{group}{Vector containing group indicators}
\item{weights}{Optional vector of sampling weights}
\item{verbose}{Logical indicating whether progress should be printed}
}

\details{
For \eqn{G} groups and \eqn{I} items, item loadings \eqn{\lambda_{ig0}}
and intercepts \eqn{\nu_{ig0}} are available and have been estimated
in a 1-dimensional factor analysis assuming a standardized factor.

The alignment procedure searches means \eqn{\alpha_{g0}}
and standard deviations \eqn{\psi_{g0}} using an alignment
optimization function \eqn{F}. This function is defined as
\deqn{F=\sum_i \sum_{ g_1 < g_2} w_{i,g1} w_{i,g2}
    f_\lambda( \lambda_{i g_1,1} - \lambda_{i g_2,1} )
    + \sum_i \sum_{ g_1 < g_2} w_{i,g1} w_{i,g2}
f_\nu( \nu_{i g_1,1} - \nu_{i g_2,1} ) }
where the aligned item parameters \eqn{\lambda_{i g,1}}
and \eqn{\nu_{i g,1}} are defined such that

\deqn{ \lambda_{i g,1}=\lambda_{i g 0} / \psi_{g0}
    \qquad \mbox{and} \qquad
    \nu_{i g,1}=\nu_{i g 0} -  \alpha_{g0} \lambda_{ig0} / \psi_{g0}
            }
and the optimization functions are defined as
\deqn{ f_\lambda (x)=| x/ a_\lambda | ^{p_\lambda}
\approx [ ( x/ a_\lambda )^2 + \varepsilon ]^{p_\lambda / 2}
    \qquad \mbox{and} \qquad
    f_\nu (x)=|  x/ a_\nu ]^{p_\nu}
    \approx [ ( x/ a_\nu )^2 + \varepsilon ]^{p_\nu / 2}
            }
using a small \eqn{ \varepsilon > 0} (e.g. .001) to obtain
a differentiable optimization function. For \eqn{p_\nu=0} or \eqn{p_\lambda=0}, the
optimization function essentially counts the number of different parameter
and mimicks a \eqn{L_0} penalty which is zero iff the argument is zero
and one otherwise. It is approximated by
\deqn{f(x)=x^2 (x^2 + \varepsilon )^{-1} }
(O'Neill & Burke, 2023).

For identification reasons, the product \eqn{\Pi_g \psi_{g0}} (\code{meth}=0,0.5)
of all group standard deviations or \eqn{\psi_1} (\code{meth}=1,2)
is set to one. The mean
\eqn{\alpha_{g0}} of the first group is set to zero (\code{meth}=0.5,1,2) or
a penalty function is added to the linking function (\code{meth}=0).

Note that Asparouhov and Muthen (2014) use \eqn{a_\lambda=a_\nu=1}
(which can be modified in \code{align.scale})
and \eqn{p_\lambda=p_\nu=0.5} (which can be modified in \code{align.pow}).
In case of \eqn{p_\lambda=2}, the penalty is approximately
\eqn{f_\lambda(x)=x^2 }, in case of \eqn{p_\lambda=0.5}
it is approximately \eqn{f_\lambda(x)=\sqrt{|x|} }. Note that \pkg{sirt} used a
different parametrization in versions up to 3.5. The \eqn{p} parameters have to be halved
for consistency with previous versions (e.g., the Asparouhov & Muthen parametrization
corresponds to \eqn{p=.25}; see also Fischer & Karl, 2019, for an application of
the previous parametrization).

Effect sizes of approximate invariance based on \eqn{R^2} have
been proposed by Asparouhov and Muthen (2014). These are
calculated separately for item loading and intercepts, resulting
in \eqn{R^2_\lambda} and \eqn{R^2_\nu} measures which are
included in the output \code{es.invariance}. In addition,
the average correlation of aligned item parameters among groups (\code{rbar})
is reported.

\emph{Metric invariance} means that all aligned item loadings \eqn{\lambda_{ig,1}}
are equal across groups and therefore \eqn{R^2_\lambda=1}.
\emph{Scalar invariance} means that all aligned item loadings
\eqn{\lambda_{ig,1}} and aligned item intercepts \eqn{\nu_{ig,1}} are
equal across groups and therefore \eqn{R^2_\lambda=1} and \eqn{R^2_\nu=1}
(see Vandenberg & Lance, 2000).
}

\value{
A list with following entries
\item{pars}{Aligned distribution parameters}
\item{itempars.aligned}{Aligned item parameters for all groups}
\item{es.invariance}{Effect sizes of approximate invariance}
\item{lambda.aligned}{Aligned \eqn{ \lambda_{i g,1}} parameters}
\item{lambda.resid}{Residuals of \eqn{ \lambda_{i g,1}} parameters}
\item{nu.aligned}{Aligned \eqn{ \nu_{i g,1}} parameters}
\item{nu.resid}{Residuals of \eqn{ \nu_{i g,1}} parameters}
\item{Niter}{Number of iterations for \eqn{f_\lambda} and
\eqn{f_\nu} optimization functions}
\item{fopt}{Minimum optimization value}
\item{align.scale}{Used alignment scale parameters}
\item{align.pow}{Used alignment power parameters}
\item{vcov}{Estimated variance matrix of aligned means and
standard deviations}
\item{\dots}{More values}
}


\references{
Asparouhov, T., & Muthen, B. (2014). Multiple-group factor analysis alignment.
\emph{Structural Equation Modeling, 21}(4), 1-14.
\doi{10.1080/10705511.2014.919210}

Byrne, B. M., & van de Vijver, F. J. R. (2017). The maximum likelihood alignment
approach to testing for approximate measurement invariance:
A paradigmatic cross-cultural application. \emph{Psicothema, 29}(4), 539-551.
\doi{10.7334/psicothema2017.178}

DeMars, C. E. (2020). Alignment as an alternative to anchor purification in DIF analyses.
\emph{Structural Equation Modeling, 27}(1), 56-72.
\doi{10.1080/10705511.2019.1617151}

Finch, W. H. (2016). Detection of differential item functioning for more
than two groups: A Monte Carlo comparison of methods.
\emph{Applied Measurement in Education, 29},(1), 30-45,
\doi{10.1080/08957347.2015.1102916}

Fischer, R., & Karl, J. A. (2019). A primer to (cross-cultural) multi-group invariance
testing possibilities in R.
\emph{Frontiers in Psychology | Cultural Psychology, 10}:1507.
\doi{10.3389/fpsyg.2019.01507}

Flake, J. K., & McCoach, D. B. (2018). An investigation of the alignment method with
polytomous indicators under conditions of partial measurement invariance.
\emph{Structural Equation Modeling, 25}(1), 56-70.
\doi{10.1080/10705511.2017.1374187}

Kim, E. S., Cao, C., Wang, Y., & Nguyen, D. T. (2017). Measurement invariance testing
with many groups: A comparison of five approaches.
\emph{Structural Equation Modeling, 24}(4), 524-544.
\doi{10.1080/10705511.2017.1304822}

Marsh, H. W., Guo, J., Parker, P. D., Nagengast, B., Asparouhov, T., Muthen, B.,
& Dicke, T. (2018). What to do when scalar invariance fails: The extended alignment
method for multi-group factor analysis comparison of latent means across many groups.
\emph{Psychological Methods, 23}(3), 524-545.
doi: 10.1037/met0000113

Muthen, B., & Asparouhov, T. (2014). IRT studies of many groups: The alignment method.
\emph{Frontiers in Psychology | Quantitative Psychology and Measurement, 5}:978.
\doi{10.3389/fpsyg.2014.00978}

Muthen, B., & Asparouhov, T. (2018). Recent methods for the study of measurement
invariance with many groups: Alignment and random effects.
\emph{Sociological Methods & Research, 47}(4), 637-664.
\doi{10.1177/0049124117701488}

O'Neill, M., & Burke, K. (2023). Variable selection using a smooth information criterion
for distributional regression models. \emph{Statistics and Computing, 33}(3), 71.
\doi{10.1007/s11222-023-10204-8}

Pokropek, A., Davidov, E., & Schmidt, P. (2019). A Monte Carlo simulation study to
assess the appropriateness of traditional and newer approaches to test for
measurement invariance. \emph{Structural Equation Modeling, 26}(5), 724-744.
\doi{10.1080/10705511.2018.1561293}

Vandenberg, R. J., & Lance, C. E. (2000). A review and synthesis of the
measurement invariance literature: Suggestions, practices, and
recommendations for organizational research. \emph{Organizational Research
Methods, 3}, 4-70.
\doi{10.1177/109442810031002}s
}

%Oelker, M. R., & Tutz, G. (2017). A uniform framework for the combination of penalties in
%generalized structured models. \emph{Advances in Data Analysis and Classification, 11}(1),
%97-120.
%\doi{10.1007/s11634-015-0205-y}

%\author{
%Alexander Robitzsch
%}

%\note{
%vvv
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For IRT linking see also \code{\link{linking.haberman}} or
\code{\link[TAM:tam.linking]{TAM::tam.linking}}.

For modeling random item effects for loadings and intercepts
see \code{\link{mcmc.2pno.ml}}.
}

\examples{
#############################################################################
# EXAMPLE 1: Item parameters cultural activities
#############################################################################

data(data.activity.itempars, package="sirt")
lambda <- data.activity.itempars$lambda
nu <- data.activity.itempars$nu
Ng <-  data.activity.itempars$N
wgt <- matrix( sqrt(Ng), length(Ng), ncol(nu) )

#***
# Model 1: Alignment using a quadratic loss function
mod1 <- sirt::invariance.alignment( lambda, nu, wgt, align.pow=c(2,2) )
summary(mod1)

#****
# Model 2: Different powers for alignment
mod2 <- sirt::invariance.alignment( lambda, nu, wgt,  align.pow=c(.5,1),
              align.scale=c(.95,.95))
summary(mod2)

# compare means from Models 1 and 2
plot( mod1$pars$alpha0, mod2$pars$alpha0, pch=16,
    xlab="M (Model 1)", ylab="M (Model 2)", xlim=c(-.3,.3), ylim=c(-.3,.3) )
lines( c(-1,1), c(-1,1), col="gray")
round( cbind( mod1$pars$alpha0, mod2$pars$alpha0 ), 3 )
round( mod1$nu.resid, 3)
round( mod2$nu.resid,3 )

# L0 penalty
mod2b <- sirt::invariance.alignment( lambda, nu, wgt,  align.pow=c(0,0),
              align.scale=c(.3,.3))
summary(mod2b)

#****
# Model 3: Low powers for alignment of scale and power
# Note that setting increment.factor larger than 1 seems necessary
mod3 <- sirt::invariance.alignment( lambda, nu, wgt, align.pow=c(.5,.75),
            align.scale=c(.55,.55), psi0.init=mod1$psi0, alpha0.init=mod1$alpha0 )
summary(mod3)

# compare mean and SD estimates of Models 1 and 3
plot( mod1$pars$alpha0, mod3$pars$alpha0, pch=16)
plot( mod1$pars$psi0, mod3$pars$psi0, pch=16)

# compare residuals for Models 1 and 3
# plot lambda
plot( abs(as.vector(mod1$lambda.resid)), abs(as.vector(mod3$lambda.resid)),
      pch=16, xlab="Residuals lambda (Model 1)",
      ylab="Residuals lambda (Model 3)", xlim=c(0,.1), ylim=c(0,.1))
lines( c(-3,3),c(-3,3), col="gray")
# plot nu
plot( abs(as.vector(mod1$nu.resid)), abs(as.vector(mod3$nu.resid)),
      pch=16, xlab="Residuals nu (Model 1)", ylab="Residuals nu (Model 3)",
      xlim=c(0,.4),ylim=c(0,.4))
lines( c(-3,3),c(-3,3), col="gray")

\dontrun{
#############################################################################
# EXAMPLE 2: Comparison 4 groups | data.inv4gr
#############################################################################

data(data.inv4gr)
dat <- data.inv4gr
miceadds::library_install("semTools")

model1 <- "
    F=~ I01 + I02 + I03 + I04 + I05 + I06 + I07 + I08 + I09 + I10 + I11
    F ~~ 1*F
    "

res <- semTools::measurementInvariance(model1, std.lv=TRUE, data=dat, group="group")
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
lambda <- matrix( ipars[ ipars$op=="=~", "est"], nrow=4,  byrow=TRUE)
colnames(lambda) <- colnames(dat)[-1]
# extract nu's
nu <- matrix( ipars[ ipars$op=="~1"  & ipars$se !=0, "est" ], nrow=4,  byrow=TRUE)
colnames(nu) <- colnames(dat)[-1]

# Model 1: least squares optimization
mod1 <- sirt::invariance.alignment( lambda=lambda, nu=nu )
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
mod2 <- sirt::invariance.alignment( lambda=lambda, nu=nu, align.pow=c(.5,.5) )
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
mod1 <- sirt::invariance.alignment( lambda=lambda, nu=nu, align.pow=c(2,2) )
summary(mod1)

# Model 2: sparse target function and definition of scales
mod2 <- sirt::invariance.alignment( lambda=lambda, nu=nu, control=list(trace=2) )
summary(mod2)

#############################################################################
# EXAMPLE 4: Linking with item parameters containing outliers
#############################################################################

# see Help file in linking.robust

# simulate some item difficulties in the Rasch model
I <- 38
set.seed(18785)
itempars <- data.frame("item"=paste0("I",1:I) )
itempars$study1 <- stats::rnorm( I, mean=.3, sd=1.4 )
# simulate DIF effects plus some outliers
bdif <- stats::rnorm(I, mean=.4, sd=.09) +
             (stats::runif(I)>.9 )*rep( 1*c(-1,1)+.4, each=I/2 )
itempars$study2 <- itempars$study1 + bdif
# create input for function invariance.alignment
nu <- t( itempars[,2:3] )
colnames(nu) <- itempars$item
lambda <- 1+0*nu

# linking using least squares optimization
mod1 <- sirt::invariance.alignment( lambda=lambda, nu=nu )
summary(mod1)
  ##   Group Means and Standard Deviations
  ##          alpha0 psi0
  ##   study1 -0.286    1
  ##   study2  0.286    1

# linking using powers of .5
mod2 <- sirt::invariance.alignment( lambda=lambda, nu=nu, align.pow=c(1,1) )
summary(mod2)
  ##   Group Means and Standard Deviations
  ##          alpha0 psi0
  ##   study1 -0.213    1
  ##   study2  0.213    1

# linking using powers of .25
mod3 <- sirt::invariance.alignment( lambda=lambda, nu=nu, align.pow=c(.5,.5) )
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
dat.male <- dat[ dat$female==0, substring( colnames(dat),1,1)=="M"  ]
dat.female <- dat[ dat$female==1, substring( colnames(dat),1,1)=="M"  ]

#*************************
# Model 1: Linking using the Rasch model
mod1m <- sirt::rasch.mml2( dat.male )
mod1f <- sirt::rasch.mml2( dat.female )

# create objects for invariance.alignment
nu <- rbind( mod1m$item$thresh, mod1f$item$thresh )
colnames(nu) <- mod1m$item$item
rownames(nu) <- c("male", "female")
lambda <- 1+0*nu

# mean of item difficulties
round( rowMeans(nu), 3 )

# Linking using least squares optimization
res1a <- sirt::invariance.alignment( lambda, nu, align.scale=c( .3, .5 ) )
summary(res1a)

# Linking using optimization with absolute value function (pow=.5)
res1b <- sirt::invariance.alignment( lambda, nu, align.scale=c( .3, .5 ),
                align.pow=c(1,1) )
summary(res1b)

#-- compare results with Haberman linking
I <- ncol(dat.male)
itempartable <- data.frame( "study"=rep( c("male", "female"), each=I ) )
itempartable$item <- c( paste0(mod1m$item$item),  paste0(mod1f$item$item) )
itempartable$a <- 1
itempartable$b <- c( mod1m$item$b, mod1f$item$b )
# estimate linking parameters
res1c <- sirt::linking.haberman( itempars=itempartable )

#-- results of sirt::equating.rasch
x <- itempartable[ 1:I, c("item", "b") ]
y <- itempartable[ I + 1:I, c("item", "b") ]
res1d <- sirt::equating.rasch( x, y )
round( res1d$B.est, 3 )
  ##     Mean.Mean Haebara Stocking.Lord
  ##   1     0.032   0.032         0.029

#*************************
# Model 2: Linking using the 2PL model
I <- ncol(dat.male)
mod2m <- sirt::rasch.mml2( dat.male, est.a=1:I)
mod2f <- sirt::rasch.mml2( dat.female, est.a=1:I)

# create objects for invariance.alignment
nu <- rbind( mod2m$item$thresh, mod2f$item$thresh )
colnames(nu) <- mod2m$item$item
rownames(nu) <- c("male", "female")
lambda <- rbind( mod2m$item$a, mod2f$item$a )
colnames(lambda) <- mod2m$item$item
rownames(lambda) <- c("male", "female")

res2a <- sirt::invariance.alignment( lambda, nu, align.scale=c( .3, .5 ) )
summary(res2a)

res2b <- sirt::invariance.alignment( lambda, nu, align.scale=c( .3, .5 ),
                align.pow=c(1,1) )
summary(res2b)

# compare results with Haberman linking
I <- ncol(dat.male)
itempartable <- data.frame( "study"=rep( c("male", "female"), each=I ) )
itempartable$item <- c( paste0(mod2m$item$item),  paste0(mod2f$item$item ) )
itempartable$a <- c( mod2m$item$a, mod2f$item$a )
itempartable$b <- c( mod2m$item$b, mod2f$item$b )
# estimate linking parameters
res2c <- sirt::linking.haberman( itempars=itempartable )

#############################################################################
# EXAMPLE 6: Data from Asparouhov & Muthen (2014) simulation study
#############################################################################

G <- 3  # number of groups
I <- 5  # number of items
# define lambda and nu parameters
lambda <- matrix(1, nrow=G, ncol=I)
nu <- matrix(0, nrow=G, ncol=I)

# define size of noninvariance
dif <- 1

#- 1st group: N(0,1)
lambda[1,3] <- 1+dif*.4; nu[1,5] <- dif*.5

#- 2nd group: N(0.3,1.5)
gg <- 2 ; mu <- .3; sigma <- sqrt(1.5)
lambda[gg,5] <- 1-.5*dif; nu[gg,1] <- -.5*dif
nu[gg,] <- nu[gg,] + mu*lambda[gg,]
lambda[gg,] <- lambda[gg,] * sigma

#- 3nd group: N(.8,1.2)
gg <- 3 ; mu <- .8; sigma <- sqrt(1.2)
lambda[gg,4] <- 1-.7*dif; nu[gg,2] <- -.5*dif
nu[gg,] <- nu[gg,] + mu*lambda[gg,]
lambda[gg,] <- lambda[gg,] * sigma

# define alignment scale
align.scale <- c(.2,.4)   # Asparouhov and Muthen use c(1,1)
# define alignment powers
align.pow <- c(.5,.5)   # as in Asparouhov and Muthen

#*** estimate alignment parameters
mod1 <- sirt::invariance.alignment( lambda, nu, eps=.01, optimizer="optim",
            align.scale=align.scale, align.pow=align.pow, center=FALSE )
summary(mod1)

#--- find parameter constraints for prespecified tolerance
cmod1 <- sirt::invariance_alignment_constraints(model=mod1, nu_parm_tol=.4,
            lambda_parm_tol=.2 )
summary(cmod1)

#############################################################################
# EXAMPLE 7: Similar to Example 6, but with data simulation and CFA estimation
#############################################################################

#--- data simulation

set.seed(65)
G <- 3  # number of groups
I <- 5  # number of items
# define lambda and nu parameters
lambda <- matrix(1, nrow=G, ncol=I)
nu <- matrix(0, nrow=G, ncol=I)
err_var <- matrix(1, nrow=G, ncol=I)

# define size of noninvariance
dif <- 1
#- 1st group: N(0,1)
lambda[1,3] <- 1+dif*.4; nu[1,5] <- dif*.5
#- 2nd group: N(0.3,1.5)
gg <- 2 ;
lambda[gg,5] <- 1-.5*dif; nu[gg,1] <- -.5*dif
#- 3nd group: N(.8,1.2)
gg <- 3
lambda[gg,4] <- 1-.7*dif; nu[gg,2] <- -.5*dif
#- define distributions of groups
mu <- c(0,.3,.8)
sigma <- sqrt(c(1,1.5,1.2))
N <- rep(1000,3) # sample sizes per group

#* simulate data
dat <- sirt::invariance_alignment_simulate(nu, lambda, err_var, mu, sigma, N)
head(dat)

#--- estimate CFA models
pars <- sirt::invariance_alignment_cfa_config(dat[,-1], group=dat$group)
print(pars)

#--- invariance alignment
# define alignment scale
align.scale <- c(.2,.4)
# define alignment powers
align.pow <- c(.5,.5)
mod1 <- sirt::invariance.alignment( lambda=pars$lambda, nu=pars$nu, eps=.01,
            optimizer="optim", align.scale=align.scale, align.pow=align.pow, center=FALSE)
#* find parameter constraints for prespecified tolerance
cmod1 <- sirt::invariance_alignment_constraints(model=mod1, nu_parm_tol=.4,
            lambda_parm_tol=.2 )
summary(cmod1)

#--- estimate CFA models with sampling weights

#* simulate weights
weights <- stats::runif(sum(N), 0, 2)
#* estimate models
pars2 <- sirt::invariance_alignment_cfa_config(dat[,-1], group=dat$group, weights=weights)
print(pars2$nu)
print(pars$nu)

#--- estimate one-parameter model
pars <- sirt::invariance_alignment_cfa_config(dat[,-1], group=dat$group, model="1PM")
print(pars)

#############################################################################
# EXAMPLE 8: Computation of standard errors
#############################################################################

G <- 3  # number of groups
I <- 5  # number of items
# define lambda and nu parameters
lambda <- matrix(1, nrow=G, ncol=I)
nu <- matrix(0, nrow=G, ncol=I)

# define size of noninvariance
dif <- 1

mu1 <- c(0,.3,.8)
sigma1 <- c(1,1.25,1.1)

#- 1st group
lambda[1,3] <- 1+dif*.4; nu[1,5] <- dif*.5

#- 2nd group
gg <- 2
lambda[gg,5] <- 1-.5*dif; nu[gg,1] <- -.5*dif

#- 3nd group
gg <- 3
lambda[gg,4] <- 1-.7*dif; nu[gg,2] <- -.5*dif

dat <- sirt::invariance_alignment_simulate(nu=nu, lambda=lambda, err_var=1+0*lambda,
                mu=mu1, sigma=sigma1, N=500, output="data", exact=TRUE)

#* estimate CFA
res <- sirt::invariance_alignment_cfa_config(dat=dat[,-1], group=dat$group )

#- perform invariance alignment
eps <- .001
align.pow <- 0.5*rep(1,2)
lambda <- res$lambda
nu <- res$nu
mod1 <- sirt::invariance.alignment( lambda=lambda, nu=nu, eps=eps, optimizer="optim",
             align.pow=align.pow, meth=meth, vcov=res$vcov)
# variance matrix and standard errors
mod1$vcov
sqrt(diag(mod1$vcov))

#############################################################################
# EXAMPLE 9: Comparison 2 groups for dichotomous data | data.pisaMath
#############################################################################

data(data.pisaMath)
dat <- data.pisaMath$data
library("lavaan")

model1 <- "
    F=~ M192Q01 + M406Q01 + M406Q02 + M423Q01 + M496Q01 + M496Q02 + M564Q01 +
         M564Q02 + M571Q01 + M603Q01 + M603Q02
    "

fit.configural <- lavaan::cfa(model1, data=dat, group="female",
                  ordered=TRUE, std.lv=TRUE, parameterization="theta")
lavaan::summary(fit.configural, standardized=TRUE)

# extract item parameters separate group analyses
ipars <- lavaan::parameterEstimates(fit.configural)
# extract lambda's: groups are in rows, items in columns
lambda <- matrix( ipars[ ipars$op=="=~", "est"], nrow=2,  byrow=TRUE)
colnames(lambda) <- colnames(dat)[6:16]
# extract nu's
nu <- matrix( ipars[ ipars$op=="|" & ipars$se !=0, "est" ], nrow=2,  byrow=TRUE)
colnames(nu) <- colnames(dat)[6:16]

# Model 1: apply invariance alignment
mod1 <- sirt::invariance.alignment( lambda=lambda, nu=nu )
summary(mod1)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Linking}
%% \keyword{Approximate Invariance}
%% \keyword{Alignment}
%% \keyword{summary}

