%% File Name: linking.haebara.Rd
%% File Version: 0.17

\name{linking.haebara}
\alias{linking.haebara}
\alias{summary.linking.haebara}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Haebara Linking of the 2PL Model for Multiple Studies
}


\description{
The function \code{linking.haebara} is a generalization of Haebara linking
of the 2PL model to multiple groups (or multiple studies; see Battauz, 2017,
for a similar approach). The optimization estimates transformation parameters for
means and standard deviations of the groups and joint item parameters.
The function allows tow different distance functions \code{dist="L2"} and
\code{dist="L1"} where the latter is a robustified version of
Haebara linking (see Details; He, Cui, & Osterlind, 2015;
Hu, Rogers, & Vukmirovic, 2008).
}


\usage{
linking.haebara(itempars, dist="L2", theta=seq(-4,4, length=61),
        optimizer="optim", center=FALSE, eps=1e-5, use_rcpp=TRUE, ...)

\method{summary}{linking.haebara}(object, digits=3, file=NULL, ...)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
\item{itempars}{
A data frame with four or five columns. The first four columns contain
in the order: study name, item name, \eqn{a} parameter, \eqn{b} parameter.
The fifth column is an optional weight for every item and every study. See
\code{\link{linking.haberman}} for a function which uses the same argument.
}
\item{dist}{Distance function. Options are \code{"L2"} for squared loss and
\code{"L1"} for absolute value loss.}
\item{theta}{Grid of theta points for 2PL item response functions}
\item{optimizer}{Name of the optimizer chosen for alignment. Options are
   \code{"optim"} (using \code{\link[stats:optim]{stats::optim}})
   or \code{"nlminb"} (using \code{\link[stats:nlminb]{stats::nlminb}}).
}
\item{center}{Logical indicating whether means and standard deviations should
be centered after estimation}
\item{eps}{Small value for smooth approximation of the absolute value function}
\item{use_rcpp}{Logical indicating whether \pkg{Rcpp} is used for computation}
\item{object}{Object of class \code{linking.haabara}.}
\item{digits}{Number of digits after decimals for rounding in \code{summary}.}
\item{file}{Optional file name if \code{summary} should be sunk into a file.}
\item{\dots}{Further arguments to be passed}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETAILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



\details{
For \eqn{t=1,\ldots,T} studies, item difficulties \eqn{b_{it}} and
item slopes \eqn{a_{it}} are available. The 2PL item response functions are given by
\deqn{ logit P(X_{pi}=1| \theta_p )=a_i ( \theta_p - b_i ) }

Haebara linking compares the observed item response functions \eqn{P_{it}}
based on the equation for the logits \eqn{a_{it}(\theta - b_{it})} and the expected
item response functions \eqn{P_{it}^\ast} based on the equation for the logits
\eqn{a_i^\ast \sigma_t ( \theta - ( b_i - \mu_t)/\sigma_t )} where the joint
item parameters \eqn{a_i} and \eqn{b_i} and means \eqn{\mu_t} and standard
deviations \eqn{\sigma_t} are estimated.

Two loss functions are implemented. The quadratic loss of Haebara linking
(\code{dist="L2"}) minimizes
\deqn{f_{opt, L2}=\sum_t \sum_i \int ( P_{it} (\theta ) - P_{it}^\ast (\theta ) )^2 w(\theta)}
was originally proposed by Haebara. A robustified version (\code{dist="L1"})
uses the optimization function (He et al., 2015)
\deqn{f_{opt, L1}=\sum_t \sum_i \int | P_{it} (\theta ) - P_{it}^\ast (\theta ) |  w(\theta)}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\value{
A list with following entries
  \item{pars}{Estimated means and standard deviations (transformation parameters)}
  \item{item}{Estimated joint item parameters}
  \item{a.orig}{Original \eqn{a_{it}} parameters}
  \item{b.orig}{Original \eqn{b_{it}} parameters}
  \item{a.resid}{Residual \eqn{a_{it}} parameters (DIF parameters)}
  \item{b.resid}{Residual \eqn{b_{it}} parameters (DIF parameters)}
  \item{res_optim}{Value of optimization routine}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\references{
Battauz, M. (2017). Multiple equating of separate IRT calibrations.
\emph{Psychometrika, 82}, 610-636.

He, Y., Cui, Z., & Osterlind, S. J. (2015). New robust scale transformation methods
in the presence of outlying common items. \emph{Applied Psychological Measurement, 39}(8),
613-626.

Hu, H., Rogers, W. T., & Vukmirovic, Z. (2008). Investigation of IRT-based equating
methods in the presence of outlier common items. \emph{Applied Psychological Measurement,
32}(4), 311-333.
}


%\author{
%Alexander Robitzsch
%}



%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{invariance.alignment}} and \code{\link{linking.haberman}}
for alternative linking methods in the \pkg{sirt} package. See also
\code{TAM::tam.linking} in the \pkg{TAM} package for more linking functionality
and the \R packages \pkg{plink}, \pkg{equateIRT}, \pkg{equateMultiple} and
\pkg{SNSequate}.
}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Robust linking methods in the presence of outliers
#############################################################################

#** simulate data
I <- 10
a <- seq(.9, 1.1, len=I)
b <- seq(-2, 2, len=I)

#- define item parameters
item_names <- paste0("I",100+1:I)
# th=SIG*TH+MU=> logit(p)=a*(SIG*TH+MU-b)=a*SIG*(TH-(-MU)/SIG-b/SIG)
d1 <- data.frame( study="S1", item=item_names, a=a, b=b )
mu <- .5; sigma <- 1.3
d2 <- data.frame( study="S2", item=item_names, a=a*sigma, b=(b-mu)/sigma )
mu <- -.3; sigma <- .7
d3 <- data.frame( study="S3", item=item_names, a=a*sigma, b=(b-mu)/sigma )

#- define DIF effect
# dif <- 0  # no DIF effects
dif <- 1
d2[4,"a"] <- d2[4,"a"] * (1-.8*dif)
d3[5,"b"] <- d3[5,"b"] - 2*dif
itempars <- rbind(d1, d2, d3)

#* Haebara linking non-robust
mod1 <- sirt::linking.haebara( itempars, dist="L2", control=list(trace=2) )
summary(mod1)

#* Haebara linking robust
mod2 <- sirt::linking.haebara( itempars, dist="L1", control=list(trace=2) )
summary(mod2)

#* Haberman linking non-robust
mod3 <- sirt::linking.haberman(itempars)
summary(mod3)

#* Haberman linking robust
mod4 <- sirt::linking.haberman(itempars, a_trim=.25, b_trim=.5)
summary(mod4)

#* compare transformation parameters (means and standard deviations)
mod1$pars
mod2$pars
mod3$transf.personpars
mod4$transf.personpars
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{Linking}
%\keyword{Equating}
