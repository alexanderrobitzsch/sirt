%% File Name: equating.rasch.Rd
%% File Version: 1.27

\name{equating.rasch}
\alias{equating.rasch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Equating in the Generalized Logistic Rasch Model
}
\description{
This function does the linking in the generalized
logistic item response model. Only item difficulties (\eqn{b}
item parameters) are allowed. Mean-mean linking and the methods 
of Haebara and Stocking-Lord are implemented (Kolen & Brennan, 2004).
}
\usage{
equating.rasch(x, y, theta = seq(-4, 4, len = 100), 
       alpha1 = 0, alpha2 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Matrix with two columns: First column items, second column item
difficulties
}
  \item{y}{
Matrix with two columns: First columns item, second column item
difficulties
}
  \item{theta}{
Vector of theta values at which the linking functions
should be evaluated. If a weighting according to a prespecified normal 
distribution \eqn{N( \mu,\sigma^2)} is aimed, then choose
\code{theta=stats::qnorm( seq(.001 , .999 , len=100) , mean=mu, sd=sigma)}
}
  \item{alpha1}{
Fixed \eqn{\alpha_1} parameter in the generalized item
response model
}
  \item{alpha2}{
Fixed \eqn{\alpha_2} parameter in the generalized item
response model
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{B.est}{Estimated linking constants according to the methods 
\code{Mean.Mean} (Mean-mean linking), \code{Haebara} (Haebara method)
and \code{Stocking.Lord} (Stocking-Lord method).}
\item{descriptives}{Descriptives of the linking. The linking error 
(\code{linkerror}) is calculated under the assumption of simple 
random sampling of items}
\item{anchor}{Original and transformed item parameters of anchor items}
\item{transf.par}{Original and transformed item parameters of all items}
}
\references{
Kolen, M. J., & Brennan, R. L. (2004). \emph{Test Equating, Scaling, and Linking:
Methods and Practices}. New York: Springer.
}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For estimating standard errors (due to inference with respect to
the item domain) of this procedure see \code{\link{equating.rasch.jackknife}}.

For linking several studies see \code{\link{linking.haberman}} or 
\code{\link{invariance.alignment}}.

A robust alternative to mean-mean linking is implemented in
\code{\link{linking.robust}}.

For linking under more general item response models
see the \pkg{plink} package.
}


\examples{
#############################################################################
# EXAMPLE 1: Linking item parameters of the PISA study
#############################################################################

data(data.pisaPars)
pars <- data.pisaPars

# linking the two studies with the Rasch model
mod <- sirt::equating.rasch(x=pars[,c("item","study1")], y=pars[,c("item","study2")])
  ##   Mean.Mean    Haebara Stocking.Lord
  ## 1   0.08828 0.08896269    0.09292838

\dontrun{
#*** linking using the plink package
# The plink package is not available on CRAN anymore.
# You can download the package with
# utils::install.packages("plink", repos = "http://www2.uaem.mx/r-mirror")
library(plink)    
I <- nrow(pars)
pm <- plink::as.poly.mod(I)
# linking parameters
plink.pars1 <- list( "study1" = data.frame( 1 , pars$study1 , 0 ) ,
                     "study2" = data.frame( 1 , pars$study2 , 0 ) )
      # the parameters are arranged in the columns:
      # Discrimination, Difficulty, Guessing Parameter
# common items
common.items <- cbind("study1"=1:I,"study2"=1:I)
# number of categories per item
cats.item <- list( "study1"=rep(2,I), "study2"=rep(2,I))
# convert into plink object
x <- plink::as.irt.pars( plink.pars1, common.items , cat= cats.item, 
          poly.mod=list(pm,pm))
# linking using plink: first group is reference group
out <- plink::plink(x, rescale="MS", base.grp=1, D=1.7)
# summary for linking
summary(out)
  ##   -------  group2/group1*  -------
  ##   Linking Constants
  ##   
  ##                        A         B
  ##   Mean/Mean     1.000000 -0.088280
  ##   Mean/Sigma    1.000000 -0.088280
  ##   Haebara       1.000000 -0.088515
  ##   Stocking-Lord 1.000000 -0.096610
# extract linked parameters
pars.out <- plink::link.pars(out)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Equating}
\keyword{Linking}
