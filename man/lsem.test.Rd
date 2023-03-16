%% File Name: lsem.test.Rd
%% File Version: 0.07

\name{lsem.test}
\alias{lsem.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test a Local Structural Equation Model Based on Bootstrap
}
\description{
Performs global and parameter tests for a fitted local structural equation model.
The LSEM must have been fitted and bootstrap estimates of the LSEM model must be
available for statistical inference. The hypothesis of a constant parameter is tested
by means of a Wald test. Moreover, regression functions can be specified and tested
if these are specified in the argument \code{models}.
}
\usage{
lsem.test(mod, bmod, models=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{
Fitted LSEM object
}
  \item{bmod}{
Fitted LSEM bootstrap object. The argument \code{bmod} can also be missing.
}
  \item{models}{
List of model formulas for named LSEM model parameters
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
List with following entries

\item{wald_test_global}{Global Wald test for model parameters}
\item{test_models}{Output for fitted regression models}
\item{parameters}{Original model parameters after fitting (i.e., smoothing) a particular
parameter using a regression model specified in \code{models}.}
\item{parameters_boot}{Bootstrapped model parameters after fitting (i.e., smoothing)
a particular parameter using a regression model specified in \code{models}.}
}

%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{lsem.estimate}}.
}



\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: data.lsem01 | Age differentiation and tested models
#############################################################################

data(data.lsem01, package="sirt")
dat <- data.lsem01

# specify lavaan model
lavmodel <- "
        F=~ v1+v2+v3+v4+v5
        F ~~ 1*F
    "

# define grid of moderator variable age
moderator.grid <- seq(4,23,1)

#-- estimate LSEM with bandwidth 2
mod <- sirt::lsem.estimate( dat, moderator="age", moderator.grid=moderator.grid,
               lavmodel=lavmodel, h=2, std.lv=TRUE)
summary(mod1)

#-- bootstrap model
bmod <- sirt::lsem.bootstrap(mod, R=200)

#-- test models
models <- list( "F=~v1"=y ~ m + I(m^2),
                "F=~v2"=y ~ I( splines::bs(m, df=4) ) )
tmod <- sirt::lsem.test(mod=mod, bmod=bmod, models=models)
str(tmod)
sirt::print_digits(wald_test_global, 3)
sirt::print_digits(test_models, 3)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
