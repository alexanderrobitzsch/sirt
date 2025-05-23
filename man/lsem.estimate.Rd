%% File Name: lsem.estimate.Rd
%% File Version: 0.777

\name{lsem.estimate}
\alias{lsem.estimate}
\alias{summary.lsem}
\alias{plot.lsem}
\alias{lsem.MGM.stepfunctions}
\alias{lsem_local_weights}
\alias{lsem.bootstrap}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Local Structural Equation Models (LSEM)
}
\description{
Local structural equation models (LSEM) are structural equation models (SEM)
which are evaluated for each value of a pre-defined moderator variable
(Hildebrandt et al., 2009, 2016).
As in nonparametric regression models, observations near a focal point - at
which the model is evaluated - obtain higher weights, far distant observations
obtain lower weights. The LSEM can be specified by making use of \pkg{lavaan} syntax.
It is also possible to specify a discretized version of LSEM in which
values of the moderator are grouped and a multiple group SEM is specified.
The LSEM can be tested by employing a permutation test, see
\code{\link{lsem.permutationTest}}.

The function \code{lsem.MGM.stepfunctions} outputs stepwise functions
for a multiple group model evaluated at a grid of focal points of the
moderator, specified in \code{moderator.grid}.

The argument \code{pseudo_weights} provides an ad hoc solution to estimate
an LSEM for any model which can be fitted in \pkg{lavaan}.

It is also possible to constrain some of the parameters along the values
of the moderator in a joint estimation approach (\code{est_joint=TRUE}). Parameter
names can be specified which are assumed to be invariant (in \code{par_invariant}).
In addition, linear or quadratic constraints can be imposed on
parameters (\code{par_linear} or \code{par_quadratic}).

Statistical inference in case of joint estimation (but also for separate estimation)
can be conducted via bootstrap using the function \code{lsem.bootstrap}.
Bootstrap at the level of a cluster identifier is allowed (argument \code{cluster}).
}
\usage{
lsem.estimate(data, moderator, moderator.grid, lavmodel, type="LSEM", h=1.1, bw=NULL,
    residualize=TRUE, fit_measures=c("rmsea", "cfi", "tli", "gfi", "srmr"),
    standardized=FALSE, standardized_type="std.all", lavaan_fct="sem",
    sufficient_statistics=TRUE, pseudo_weights=0,
    sampling_weights=NULL, loc_linear_smooth=TRUE, est_joint=FALSE, par_invariant=NULL,
    par_linear=NULL, par_quadratic=NULL, partable_joint=NULL, pw_linear=1,
    pw_quadratic=1, pd=TRUE, est_DIF=FALSE, se=NULL, kernel="gaussian",
    eps=1e-08, verbose=TRUE, ...)

\method{summary}{lsem}(object, file=NULL, digits=3, ...)

\method{plot}{lsem}(x, parindex=NULL, ask=TRUE, ci=TRUE, lintrend=TRUE,
       parsummary=TRUE, ylim=NULL, xlab=NULL,  ylab=NULL, main=NULL,
       digits=3, ...)

lsem.MGM.stepfunctions( object, moderator.grid )

# compute local weights
lsem_local_weights(data.mod, moderator.grid, h, sampling_weights=NULL, bw=NULL,
     kernel="gaussian", is_imputed=FALSE, Nimp=0, data=NULL, moderator=NULL)

lsem.bootstrap(object, R=100, verbose=TRUE, cluster=NULL,
     repl_design=NULL, repl_factor=NULL, use_starting_values=TRUE,
     n.core=1, cl.type="PSOCK")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data frame or a list of imputed datasets
}
  \item{moderator}{
Variable name of the moderator
}
  \item{moderator.grid}{
Focal points at which the LSEM should be evaluated. If \code{type="MGM"},
breaks are defined in this vector.
}
  \item{lavmodel}{
Specified SEM in \pkg{lavaan}.
}
\item{type}{Type of estimated model. The default is \code{type="LSEM"} which means
    that a local structural equation model is estimated.
    A multiple group model with a discretized moderator as the
    grouping variable can be estimated with \code{type="MGM"}. In this
    case, the breaks must be defined in \code{moderator.grid}.
            }
  \item{h}{
Bandwidth factor
}
\item{bw}{Optional bandwidth parameter if \code{h} should not be used}
\item{residualize}{Logical indicating whether a residualization
    should be applied. }
  \item{fit_measures}{
Vector with names of fit measures following the labels in \pkg{lavaan}
}
\item{standardized}{Optional logical indicating whether
 standardized solution should be included as parameters in
 the output using the
\code{\link[lavaan:standardizedSolution]{lavaan::standardizedSolution}}
function. Standardized parameters are labeled as \code{std__}.
}
\item{standardized_type}{Type of standardization if \code{standardized=TRUE}.
The types are described in
\code{\link[lavaan:standardizedSolution]{lavaan::standardizedSolution}}.
}
\item{lavaan_fct}{String whether
\code{\link[lavaan:lavaan]{lavaan::lavaan}} (\code{lavaan_fct="lavaan"}),
\code{\link[lavaan:sem]{lavaan::sem}} (\code{lavaan_fct="sem"}),
\code{\link[lavaan:cfa]{lavaan::cfa}} (\code{lavaan_fct="cfa"}) or
\code{\link[lavaan:growth]{lavaan::growth}} (\code{lavaan_fct="growth"})
should be used.
}
\item{sufficient_statistics}{Logical whether sufficient statistics of weighted
means and covariances should be used for model fitting. This option
can be set to \code{sufficient_statistics=FALSE}
if the data contain missing values. Note that the option
\code{sufficient_statistics=TRUE} is
only valid for (approximate) missing completely at random (MCAR) data.
The option can only be used for continuous data.}
\item{pseudo_weights}{Integer defining a target sample size. Local weights
are multiplied by a factor which is rounded to integers.
This approach is referred as a pseudo weighting approach.
For example, using \code{pseudo_weights=30000} implies
that the sum of local weights at each focal point is \code{30000}.
}
\item{sampling_weights}{Optional vector of sampling weights}
\item{loc_linear_smooth}{Logical indicating whether local linear
smoothing should be used for computing sufficient statistics for
means and covariances. The default is \code{FALSE}.}
\item{est_joint}{Logical indicating whether LSEM should be estimated in a
joint estimation approach. This options only works wih continuous data and
sufficient statistics.}
\item{par_invariant}{Vector of invariant parameters}
\item{par_linear}{Vector of parameters with linear function}
\item{par_quadratic}{Vector of parameters with quadratic function}
\item{partable_joint}{User-defined parameter table if joint estimation is
used (\code{est_joint=TRUE}).}
\item{pw_linear}{Number of segments if piecewise linear estimation of parameters is used}
\item{pw_quadratic}{Number of segments if piecewise quadratic estimation of parameters
is used}
\item{pd}{Logical indicating whether nearest positive definite covariance matrix
should be computed if sufficient statistics are used}
\item{est_DIF}{Logical indicating whether parameters under differential
item functioning (DIF) should be additionally
computed for invariant item parameters}
\item{se}{Type of standard error used in \code{lavaan::lavaan}. If \code{NULL},
the \pkg{lavaan} default is used.}
\item{kernel}{Type of kernel function. Can be \code{"gaussian"},
\code{"uniform"} or \code{"epanechnikov"}.}
\item{is_imputed}{Logical indicating whether datasets are imputed}
\item{Nimp}{Number of imputed datasets}
\item{data}{Optional list containing imputed datasets}
\item{eps}{Minimum number for weights}
\item{verbose}{Optional logical printing information about computation progress.}
\item{object}{Object of class \code{lsem}}
\item{file}{A file name in which the summary output will be written.}
\item{digits}{Number of digits.}
\item{x}{Object of class \code{lsem}.}
\item{parindex}{Vector of indices for parameters in plot function.}
\item{ask}{A logical which asks for changing the graphic for each parameter.}
\item{ci}{Logical indicating whether confidence intervals should be plotted.}
\item{lintrend}{Logical indicating whether a linear trend should be plotted.}
\item{parsummary}{Logical indicating whether a parameter summary
  should be displayed.}
\item{ylim}{Plot parameter \code{ylim}. Can be a list, see Examples.}
\item{xlab}{Plot parameter \code{xlab}. Can be a vector.}
\item{ylab}{Plot parameter \code{ylab}. Can be a vector.}
\item{main}{Plot parameter \code{main}. Can be a vector.}
  \item{\dots}{
Further arguments to be passed to \code{\link[lavaan:sem]{lavaan::sem}} or
\code{\link[lavaan:lavaan]{lavaan::lavaan}}.
}
\item{data.mod}{Observed values of the moderator}
\item{R}{Number of bootstrap samples}
\item{cluster}{Optional variable name for bootstrap at the level of a cluster identifier}
\item{repl_design}{Optional matrix containing replication weights for computation of
standard errors. Note that sampling weights have to be already included in
\code{repl_design}.}
\item{repl_factor}{Replication factor in variance formula for statistical
inference, e.g., 0.05 in PISA.}
\item{use_starting_values}{Logical indicating whether starting values should be
used from the original sample}
\item{n.core}{A scalar indicating the number of cores that should be used.}
\item{cl.type}{The cluster type.
Default value is \code{"PSOCK"}. Posix machines (Linux, Mac) generally benefit
from much faster cluster computation if type is set to \code{type="FORK"}.
}
}

%\details{
%EXPLAIN
%}

\value{
List with following entries
\item{parameters}{Data frame with all parameters estimated at focal points of
    moderator. Bias-corrected estimates under boostrap can be found in
    the column \code{est_bc}.}
\item{weights}{Data frame with weights at each focal point}
\item{parameters_summary}{Summary table for estimated parameters}
\item{parametersM}{Estimated parameters in matrix form. Parameters are in
columns and values of the grid of the moderator are in rows.}
\item{bw}{Used bandwidth}
\item{h}{Used bandwidth factor}
\item{N}{Sample size}
\item{moderator.density}{Estimated frequencies and effective sample size for
    moderator at focal points}
\item{moderator.stat}{Descriptive statistics for moderator}
\item{moderator}{Variable name of moderator}
\item{moderator.grid}{Used grid of focal points for moderator}
\item{moderator.grouped}{Data frame with informations about grouping of
moderator if \code{type="MGM"}.}
\item{residualized.intercepts}{Estimated intercept functions used for
residualization.}
\item{lavmodel}{Used lavaan model}
\item{data}{Used data frame, possibly residualized if \code{residualize=TRUE}}
\item{model_parameters}{Model parameters in LSEM}
\item{parameters_boot}{Parameter values in each bootstrap sample
(for \code{lsem.bootstrap})}
\item{fitstats_joint_boot}{Fit statistics in each bootstrap sample
(for \code{lsem.bootstrap})}
\item{dif_effects}{Estimated item parameters under DIF}
}


\references{
Hildebrandt, A., Luedtke, O., Robitzsch, A., Sommer, C., &
Wilhelm, O. (2016). Exploring factor model parameters across continuous variables
with local structural equation models.
\emph{Multivariate Behavioral Research, 51}(2-3), 257-278.
\doi{10.1080/00273171.2016.1142856}

Hildebrandt, A., Wilhelm, O., & Robitzsch, A. (2009). Complementary and
competing factor analytic approaches for the investigation of measurement invariance.
\emph{Review of Psychology, 16}, 87-102.
}


\author{
Alexander Robitzsch, Oliver Luedtke, Andrea Hildebrandt
}

% \note{
% }


%% Note that the specified model in \code{lavaan} syntax should be explicitly saved as
%% the object \code{lavmodel} in the global envirnoment as it is done in the examples.


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{lsem.permutationTest}} for conducting a permutation test
and \code{\link{lsem.test}} for applying a Wald test to a bootstrapped LSEM model.
}

\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: data.lsem01 | Age differentiation
#############################################################################

data(data.lsem01, package="sirt")
dat <- data.lsem01

# specify lavaan model
lavmodel <- "
        F=~ v1+v2+v3+v4+v5
        F ~~ 1*F"

# define grid of moderator variable age
moderator.grid <- seq(4,23,1)

#********************************
#*** Model 1: estimate LSEM with bandwidth 2
mod1 <- sirt::lsem.estimate( dat, moderator="age", moderator.grid=moderator.grid,
               lavmodel=lavmodel, h=2, std.lv=TRUE)
summary(mod1)
plot(mod1, parindex=1:5)

# perform permutation test for Model 1
pmod1 <- sirt::lsem.permutationTest( mod1, B=10 )
          # only for illustrative purposes the number of permutations B is set
          # to a low number of 10
summary(pmod1)
plot(pmod1, type="global")

#* perform permutation test with parallel computation
pmod1a <- sirt::lsem.permutationTest( mod1, B=10, n.core=3 )
summary(pmod1a)

#** estimate Model 1 based on pseudo weights
mod1b <- sirt::lsem.estimate( dat, moderator="age", moderator.grid=moderator.grid,
               lavmodel=lavmodel, h=2, std.lv=TRUE, pseudo_weights=50 )
summary(mod1b)

#** estimation with sampling weights

# generate random sampling weights
set.seed(987)
weights <- stats::runif(nrow(dat), min=.4, max=3 )
mod1c <- sirt::lsem.estimate( dat, moderator="age", moderator.grid=moderator.grid,
               lavmodel=lavmodel, h=2, sampling_weights=weights)
summary(mod1c)

#********************************
#*** Model 2: estimate multiple group model with 4 age groups

# define breaks for age groups
moderator.grid <- seq( 3.5, 23.5, len=5) # 4 groups
# estimate model
mod2 <- sirt::lsem.estimate( dat, moderator="age", moderator.grid=moderator.grid,
           lavmodel=lavmodel, type="MGM", std.lv=TRUE)
summary(mod2)

# output step functions
smod2 <- sirt::lsem.MGM.stepfunctions( object=mod2, moderator.grid=seq(4,23,1) )
str(smod2)

#********************************
#*** Model 3: define standardized loadings as derived variables

# specify lavaan model
lavmodel <- "
        F=~ a1*v1+a2*v2+a3*v3+a4*v4
        v1 ~~ s1*v1
        v2 ~~ s2*v2
        v3 ~~ s3*v3
        v4 ~~ s4*v4
        F ~~ 1*F
        # standardized loadings
        l1 :=a1 / sqrt(a1^2 + s1 )
        l2 :=a2 / sqrt(a2^2 + s2 )
        l3 :=a3 / sqrt(a3^2 + s3 )
        l4 :=a4 / sqrt(a4^2 + s4 )
        "
# estimate model
mod3 <- sirt::lsem.estimate( dat, moderator="age", moderator.grid=moderator.grid,
               lavmodel=lavmodel, h=2, std.lv=TRUE)
summary(mod3)
plot(mod3)

#********************************
#*** Model 4: estimate LSEM and automatically include standardized solutions

lavmodel <- "
        F=~ 1*v1+v2+v3+v4
        F ~~ F"
mod4 <- sirt::lsem.estimate( dat, moderator="age", moderator.grid=moderator.grid,
               lavmodel=lavmodel, h=2, standardized=TRUE)
summary(mod4)
# permutation test (use only few permutations for testing purposes)
pmod1 <- sirt::lsem.permutationTest( mod4, B=3 )

#**** compute LSEM local weights
wgt <- sirt::lsem_local_weights(data.mod=dat$age, moderator.grid=moderator.grid,
             h=2)$weights
print(str(weights))

#********************************
#*** Model 5: invariance parameter constraints and other constraints

lavmodel <- "
        F=~ 1*v1+v2+v3+v4
        F ~~ F"
moderator.grid <- seq(4,23,4)

#- estimate model without constraints
mod5a <- sirt::lsem.estimate( dat, moderator="age", moderator.grid=moderator.grid,
               lavmodel=lavmodel, h=2, standardized=TRUE)
summary(mod5a)
# extract parameter names
mod5a$model_parameters

#- invariance constraints on residual variances
par_invariant <- c("F=~v2","v2~~v2")
mod5b <- sirt::lsem.estimate( dat, moderator="age", moderator.grid=moderator.grid,
               lavmodel=lavmodel, h=2, standardized=TRUE, par_invariant=par_invariant)
summary(mod5b)

#- bootstrap for statistical inference
bmod5b <- sirt::lsem.bootstrap(mod5b, R=100)
# inspect parameter values and standard errors
bmod5b$parameters

#- bootstrap using parallel computing (i.e., multiple cores)
bmod5ba <- sirt::lsem.bootstrap(mod5b, R=100, n.core=3)

#- user-defined replication design
R <- 100    # bootstrap samples
N <- nrow(dat)
repl_design <- matrix(0, nrow=N, ncol=R)
for (rr in 1:R){
    indices <- sort( sample(1:N, replace=TRUE) )
    repl_design[,rr] <- sapply(1:N, FUN=function(ii){ sum(indices==ii) } )
}
head(repl_design)
bmod5b1 <- sirt::lsem.bootstrap(mod5a, repl_design=repl_design, repl_factor=1/R)

#- compare model mod5b with joint estimation without constraints
mod5c <- sirt::lsem.estimate( dat, moderator="age", moderator.grid=moderator.grid,
               lavmodel=lavmodel, h=2, standardized=TRUE, est_joint=TRUE)
summary(mod5c)

#- linear and quadratic functions
par_invariant <- c("F=~v1","v2~~v2")
par_linear <- c("v1~~v1")
par_quadratic <- c("v4~~v4")

mod5d <- sirt::lsem.estimate( dat1, moderator="age", moderator.grid=moderator.grid,
            lavmodel=lavmodel, h=2, par_invariant=par_invariant, par_linear=par_linear,
            par_quadratic=par_quadratic)
summary(mod5d)

#- user-defined constraints: step functions for parameters

# inspect parameter table (from lavaan) of fitted model
pj <- mod5d$partable_joint
#* modify parameter table for user-defined constraints
# define step function for F=~v1 which is constant on intervals 1:4 and 5:7
pj2 <- pj[ pj$con==1, ]
pj2[ c(5,6), "lhs" ] <- "p1g5"
pj2 <- pj2[ -4, ]
partable_joint <- rbind(pj1, pj2)
# estimate model with constraints
mod5e <- lsem::lsem.estimate( dat1, moderator="age", moderator.grid=moderator.grid,
             lavmodel=lavmodel, h=2, std.lv=TRUE, estimator="ML",
             partable_joint=partable_joint)
summary(mod5e)

#############################################################################
# EXAMPLE 2: data.lsem01 | FIML with missing data
#############################################################################

data(data.lsem01)
dat <- data.lsem01
# induce artifical missing values
set.seed(98)
dat[ stats::runif(nrow(dat)) < .5, c("v1")] <- NA
dat[ stats::runif(nrow(dat)) < .25, c("v2")] <- NA

# specify lavaan model
lavmodel1 <- "
        F=~ v1+v2+v3+v4+v5
        F ~~ 1*F"

# define grid of moderator variable age
moderator.grid <- seq(4,23,2)

#*** estimate LSEM with FIML
mod1 <- sirt::lsem.estimate( dat, moderator="age", moderator.grid=moderator.grid,
                lavmodel=lavmodel1, h=2, std.lv=TRUE, estimator="ML", missing="fiml")
summary(mod1)

#############################################################################
# EXAMPLE 3: data.lsem01 | WLSMV estimation
#############################################################################

data(data.lsem01)
dat <- data.lsem01

# create artificial dichotomous data
for (vv in 2:6){
dat[,vv] <- 1*(dat[,vv] > mean(dat[,vv]))
}

# specify lavaan model
lavmodel1 <- "
        F=~ v1+v2+v3+v4+v5
        F ~~ 1*F
        v1 | t1
        v2 | t1
        v3 | t1
        v4 | t1
        v5 | t1
        "

# define grid of moderator variable age
moderator.grid <- seq(4,23,2)

#*** local WLSMV estimation
mod1 <- sirt::lsem.estimate( dat, moderator="age", moderator.grid=moderator.grid,
          lavmodel=lavmodel1, h=2, std.lv=TRUE, estimator="DWLS", ordered=paste0("v",1:5),
          residualize=FALSE, pseudo_weights=10000, parameterization="THETA" )
summary(mod1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{LSEM}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
