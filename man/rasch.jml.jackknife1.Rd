%% File Name: rasch.jml.jackknife1.Rd
%% File Version: 2.12

\name{rasch.jml.jackknife1}
\alias{rasch.jml.jackknife1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Jackknifing the IRT Model Estimated by Joint Maximum Likelihood (JML)
}
\description{
Jackknife estimation is an alternative to other ad hoc proposed
methods for bias correction (Hahn & Newey, 2004).
}
\usage{
rasch.jml.jackknife1(jmlobj)
}
\arguments{
  \item{jmlobj}{
Output of \code{rasch.jml}
}
}
\details{
Note that items are used for jackknifing (Hahn & Newey, 2004). 
By default, all \eqn{I} items in the data frame are used as 
jackknife units. 
}
\value{
A list with following entries
  \item{item}{
        A data frame with item parameters
        \itemize{
            \item \code{b.JML}: {Item difficulty from JML estimation}
            \item \code{b.JMLcorr}: {Item difficulty from JML estimation by 
                    applying the correction factor \eqn{(I-1)/I}}
            \item \code{b.jack}: {Item difficulty from Jackknife estimation}
            \item \code{b.jackse}: {Standard error of Jackknife estimation 
                        for item difficulties.
                        Note that this parameter refer to the standard error 
                    with respect to item sampling}
            \item \code{b.JMLse}: {Standard error for item difficulties 
                    obtained from JML estimation}
                }
            }        
  \item{jack.itemdiff}{A matrix containing all item difficulties obtained by Jackknife}
}
\references{
Hahn, J., & Newey, W. (2004). Jackknife and analytical bias reduction for 
nonlinear panel models. \emph{Econometrica}, \bold{72}, 1295-1319.
}
\author{
Alexander Robitzsch
}

\seealso{
For JML estimation \code{\link{rasch.jml}}.

For analytical bias correction methods see \code{\link{rasch.jml.biascorr}}.
}
\examples{
#############################################################################
# EXAMPLE 1: Simulated data from the Rasch model
#############################################################################
set.seed(7655)
N <- 5000    # number of persons
I <- 11      # number of items
b <- seq( -2 , 2 , length=I )
dat <- sirt::sim.raschtype( rnorm( N ) , b )
colnames(dat) <- paste( "I" , 1:I , sep="")

# estimate the Rasch model with JML
mod <- sirt::rasch.jml( dat  )
summary(mod)

# re-estimate the Rasch model using Jackknife
mod2 <- sirt::rasch.jml.jackknife1( mod )
  ##   
  ##   Joint Maximum Likelihood Estimation 
  ##   Jackknife Estimation
  ##   11 Jackknife Units are used
  ##   |--------------------PROGRESS--------------------| 
  ##   |------------------------------------------------|
  ##   
  ##          N     p  b.JML b.JMLcorr b.jack b.jackse b.JMLse
  ##   I1  4929 0.853 -2.345    -2.131 -2.078    0.079   0.045
  ##   I2  4929 0.786 -1.749    -1.590 -1.541    0.075   0.039
  ##   I3  4929 0.723 -1.298    -1.180 -1.144    0.065   0.036
  ##   I4  4929 0.657 -0.887    -0.806 -0.782    0.059   0.035
  ##   I5  4929 0.576 -0.420    -0.382 -0.367    0.055   0.033
  ##   I6  4929 0.492  0.041     0.038  0.043    0.054   0.033
  ##   I7  4929 0.409  0.502     0.457  0.447    0.056   0.034
  ##   I8  4929 0.333  0.939     0.854  0.842    0.058   0.035
  ##   I9  4929 0.264  1.383     1.257  1.229    0.065   0.037
  ##   I10 4929 0.210  1.778     1.617  1.578    0.071   0.040
  ##   I11 4929 0.154  2.266     2.060  2.011    0.077   0.044
#-> Item parameters obtained by jackknife seem to be acceptable.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Joint maximum likelihood (JML)}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
