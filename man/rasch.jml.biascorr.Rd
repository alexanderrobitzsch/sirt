%% File Name: rasch.jml.biascorr.Rd
%% File Version: 0.17

\name{rasch.jml.biascorr}
\alias{rasch.jml.biascorr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bias Correction of Item Parameters for Joint Maximum Likelihood Estimation
in the Rasch model
}
\description{
This function computes an analytical bias correction for the Rasch
model according to the method of Arellano and Hahn (2007).
}
\usage{
rasch.jml.biascorr(jmlobj,itemfac=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{jmlobj}{
An object which is the output of the \code{\link{rasch.jml}} function
}
\item{itemfac}{Number of items which are used for bias correction.
By default it is the average number of item responses per person.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with following entries
\item{b.biascorr}{Matrix of item difficulty estimates. The column
\code{b.analytcorr1} contains item difficulties by analytical bias
correction of Method 1 in Arellano and Hahn (2007) whereas \code{b.analytcorr2}
corresponds to Method 2.}
\item{b.bias1}{Estimated bias by Method 1}
\item{b.bias2}{Estimated bias by Method 2}
\item{itemfac}{Number of items which are used as the factor for bias
    correction}
}
\references{
Arellano, M., & Hahn, J. (2007). Understanding bias in nonlinear panel
models: Some recent developments. In R. Blundell, W. Newey &
T. Persson (Eds.): \emph{Advances in Economics and
Econometrics, Ninth World Congress}, Cambridge University Press.
}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{rasch.jml.jackknife1}} for bias correction based on
Jackknife.

See also the \pkg{bife} \R package for analytical bias corrections.
}


\examples{
#############################################################################
# EXAMPLE 1: Dataset Reading
#############################################################################
data(data.read)
dat <- data( data.read )

# estimate Rasch model
mod <- sirt::rasch.jml( data.read  )

# JML with analytical bias correction
res1 <- sirt::rasch.jml.biascorr( jmlobj=mod  )
print( res1$b.biascorr , digits= 3 )
  ##        b.JML b.JMLcorr b.analytcorr1 b.analytcorr2
  ##   1  -2.0086   -1.8412        -1.908        -1.922
  ##   2  -1.1121   -1.0194        -1.078        -1.088
  ##   3  -0.0718   -0.0658        -0.150        -0.127
  ##   4   0.5457    0.5002         0.393         0.431
  ##   5  -0.9504   -0.8712        -0.937        -0.936
  ##  [...]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Joint maximum likelihood (JML)}
%\keyword{xxx}
