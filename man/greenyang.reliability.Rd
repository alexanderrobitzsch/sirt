%% File Name: greenyang.reliability.Rd
%% File Version: 1.31

\name{greenyang.reliability}
\alias{greenyang.reliability}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reliability for Dichotomous Item Response Data
Using the Method of Green and Yang (2009)
}
\description{
This function estimates the model-based reliability
of dichotomous data using the Green & Yang (2009) method.
The underlying factor model is \eqn{D}-dimensional where
the dimension \eqn{D} is specified by the argument
\code{nfactors}. The factor solution is subject to the
application of the Schmid-Leiman transformation (see Reise, 2012;
Reise, Bonifay, & Haviland, 2013; Reise, Moore, & Haviland, 2010).
}

\usage{
greenyang.reliability(object.tetra, nfactors)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object.tetra}{
Object as the output of the function \code{tetrachoric}, the
\code{fa.parallel.poly} from the \pkg{psych} package or the
\code{tetrachoric2} function (from \pkg{sirt}).
This object can also be created
as a list by the user where the tetrachoric correlation
must must be in the list entry \code{rho} and the thresholds must
be in the list entry \code{thresh}.
}
  \item{nfactors}{
Number of factors (dimensions)
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A data frame with columns:
\item{coefficient}{Name of the reliability measure. \code{omega_1} (Omega)
is the reliability estimate for the total score for dichotomous data
based on a one-factor model, \code{omega_t} (Omega Total) is the
estimate for a \eqn{D}-dimensional model. For the nested factor model,
\code{omega_h} (Omega Asymptotic) is the reliability of the general factor model,
\code{omega_ha} (Omega Hierarchical Asymptotic) eliminates item-specific
variance. The explained common variance (\code{ECV}) explained by the
common factor is based on the \eqn{D}-dimensional but does not take
item thresholds into account. The amount of explained
variance \code{ExplVar} is defined as the quotient of the first
eigenvalue of the tetrachoric correlation matrix to the
sum of all eigenvalues. The statistic \code{EigenvalRatio}
is the ratio of the first and second eigenvalue.
}
\item{dimensions}{Number of dimensions}
\item{estimate}{Reliability estimate}
}


\references{
Green, S. B., & Yang, Y. (2009). Reliability of summed item
scores using structural equation modeling: An alternative to
coefficient alpha. \emph{Psychometrika, 74}, 155-167.

Reise, S. P. (2012). The rediscovery of bifactor measurement models.
\emph{Multivariate Behavioral Research, 47}, 667-696.

Reise, S. P., Bonifay, W. E., & Haviland, M. G. (2013).
Scoring and modeling psychological measures in the presence of
multidimensionality. \emph{Journal of Personality Assessment,
95}, 129-140.

Reise, S. P., Moore, T. M., & Haviland, M. G.  (2010).
Bifactor models and rotations: Exploring the extent to which
multidimensional data yield univocal scale scores,
\emph{Journal of Personality Assessment, 92}, 544-559.
}


%\author{
%Alexander Robitzsch
%}


\note{
This function needs the \pkg{psych} package.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{f1d.irt}} for estimating the functional unidimensional
item response model.

This function uses \code{\link{reliability.nonlinearSEM}}.

See also the \code{MBESS::ci.reliability} function for estimating
reliability for polytomous item responses.
}

\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Reliability estimation of Reading dataset data.read
#############################################################################
miceadds::library_install("psych")
set.seed(789)
data( data.read )
dat <- data.read

# calculate matrix of tetrachoric correlations
dat.tetra <- psych::tetrachoric(dat)      # using tetrachoric from psych package
dat.tetra2 <- sirt::tetrachoric2(dat)       # using tetrachoric2 from sirt package

# perform parallel factor analysis
fap <- psych::fa.parallel.poly(dat, n.iter=1 )
  ##   Parallel analysis suggests that the number of factors=3
  ##   and the number of components=2

# parallel factor analysis based on tetrachoric correlation matrix
##       (tetrachoric2)
fap2 <- psych::fa.parallel(dat.tetra2$rho, n.obs=nrow(dat),  n.iter=1 )
  ## Parallel analysis suggests that the number of factors=6
  ## and the number of components=2
  ## Note that in this analysis, uncertainty with respect to thresholds is ignored.

# calculate reliability using a model with 4 factors
greenyang.reliability( object.tetra=dat.tetra, nfactors=4 )
  ##                                            coefficient dimensions estimate
  ## Omega Total (1D)                               omega_1          1    0.771
  ## Omega Total (4D)                               omega_t          4    0.844
  ## Omega Hierarchical (4D)                        omega_h          4    0.360
  ## Omega Hierarchical Asymptotic (4D)            omega_ha          4    0.427
  ## Explained Common Variance (4D)                     ECV          4    0.489
  ## Explained Variance (First Eigenvalue)          ExplVar         NA   35.145
  ## Eigenvalue Ratio (1st to 2nd Eigenvalue) EigenvalRatio         NA    2.121

# calculation of Green-Yang-Reliability based on tetrachoric correlations
#   obtained by tetrachoric2
greenyang.reliability( object.tetra=dat.tetra2, nfactors=4 )

# The same result will be obtained by using fap as the input
greenyang.reliability( object.tetra=fap, nfactors=4 ) }
}
\keyword{Reliability}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
