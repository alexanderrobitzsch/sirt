%% File Name: class.accuracy.rasch.Rd
%% File Version: 0.24

\name{class.accuracy.rasch}
\alias{class.accuracy.rasch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Classification Accuracy in the Rasch Model
}
\description{
This function computes the classification accuracy in the Rasch model
for the maximum likelihood (person parameter) estimate according
to the method of Rudner (2001).
}
\usage{
class.accuracy.rasch(cutscores, b, meantheta, sdtheta, theta.l,
    n.sims=0, seed=988)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cutscores}{
Vector of cut scores
}
  \item{b}{
Vector of item difficulties
}
  \item{meantheta}{
Mean of the trait distribution
}
  \item{sdtheta}{
Standard deviation of the trait distribution
}
  \item{theta.l}{
Discretized theta distribution
}
  \item{n.sims}{
Number of simulated persons in a data set. The default is 0
which means that no simulation is performed.
}
\item{seed}{The random seed for the simulation}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with following entries:
  \item{class.stats}{Data frame containing classification accuracy statistics. The
  column \code{agree0} refers to absolute agreement, \code{agree1} to
  the agreement of at most a difference of one level.
     }
  \item{class.prob}{Probability table of classification}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Rudner, L.M. (2001). Computing the expected proportions of misclassified examinees.
\emph{Practical Assessment, Research & Evaluation, 7}(14).
}

%\author{
%Alexander Robitzsch
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Classification accuracy of other IRT models
can be obtained with the \R package \pkg{cacIRT}.
}
\examples{
#############################################################################
# EXAMPLE 1: Reading dataset
#############################################################################
data( data.read, package="sirt")
dat <- data.read

# estimate the Rasch model
mod <- sirt::rasch.mml2( dat )

# estimate classification accuracy (3 levels)
cutscores <- c( -1, .3 )    # cut scores at theta=-1 and theta=.3
class.accuracy.rasch( cutscores=cutscores, b=mod$item$b,
           meantheta=0,  sdtheta=mod$sd.trait,
           theta.l=seq(-4,4,len=200 ), n.sims=3000)
  ##   Cut Scores
  ##   [1] -1.0  0.3
  ##
  ##   WLE reliability (by simulation)=0.671
  ##   WLE consistency (correlation between two parallel forms)=0.649
  ##
  ##   Classification accuracy and consistency
  ##              agree0 agree1 kappa consistency
  ##   analytical   0.68  0.990 0.492          NA
  ##   simulated    0.70  0.997 0.489       0.599
  ##
  ##   Probability classification table
  ##               Est_Class1 Est_Class2 Est_Class3
  ##   True_Class1      0.136      0.041      0.001
  ##   True_Class2      0.081      0.249      0.093
  ##   True_Class3      0.009      0.095      0.294
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Classification accuracy}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
