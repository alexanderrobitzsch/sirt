\name{unidim.test.csn}
\alias{unidim.test.csn}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Test for Unidimensionality of CSN 
}

\description{
This function tests whether item covariances given the sum
score are non-positive (CSN; see Junker 1993), 
i.e. for items \eqn{i} and \eqn{j} it holds that 
\deqn{ Cov( X_i , X_j | X^+ ) \le 0 }
Note that this function only works for dichotomous data.
}


\usage{
unidim.test.csn(dat, RR = 400, prop.perm = 0.75, progress = TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame with dichotomous item responses. All persons with (some)
missing responses are removed.
}
  \item{RR}{
Number of permutations used for statistical testing
}
  \item{prop.perm}{
A positive value indicating the amount of permutation in an
existing permuted data set
}
  \item{progress}{
An optional logical indicating whether computation progress
should be displayed
}
}

\details{
For each item pair \eqn{(i,j)} and a each sum score group \eqn{k}
a conditional covariance \eqn{r(i,j|k)} is calculated. Then,
the test statistic for CSN is
\deqn{ h = \sum_{k=1}^{I-1} \frac{n_k}{n} \max_{i,j} r(i,j|k) }
where \eqn{n_k} is the number of persons in score group \eqn{k}.
"'Large values"' of \eqn{h} are not in agreement with the null
hypothesis of non-positivity of conditional covariances.

The distribution of the test statistic \eqn{h} under the null
hypothesis is empirically obtained by column wise permutation
of items within all score groups. In the population, this procedure
corresponds to conditional covariances of zero. See de Gooijer and Yuan (2011)
for more details.
}

\value{
A list with following entries
\item{stat}{Value of the statistic}
\item{stat_perm}{Distribution of statistic under \eqn{H_0} of
permuted dataset}
\item{p}{The corresponding p value of the statistic}
\item{H0_quantiles}{Quantiles of the statistic under permutation
(the null hypothesis \eqn{H_0})}

}


\references{
De Gooijer, J. G., & Yuan, A. (2011). Some exact tests for
manifest properties of latent trait models.
\emph{Computational Statistics and Data Analysis}, \bold{55},
34-44.

Junker, B.W. (1993). Conditional association, essential independence, and 
monotone unidimensional item response models. 
\emph{Annals of Statistics}, \bold{21}, 1359-1378.
}
\author{
Alexander Robitzsch
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%MH statistic conditional association
%}


\examples{
#############################################################################
# EXAMPLE 1: Dataset data.read
#############################################################################
	
data(data.read)
dat <- data.read
set.seed(778)
res <- sirt::unidim.test.csn( dat )
  ##  CSN Statistic = 0.04737 , p = 0.02 

\dontrun{
#############################################################################
# EXAMPLE 2: CSN statistic for two-dimensional simulated data
#############################################################################
	
set.seed(775)
N <- 2000 
I <- 30   # number of items
rho <- .60   # correlation between 2 dimensions
t0 <- stats::rnorm(N)
t1 <- sqrt(rho)*t0 + sqrt(1-rho)*stats::rnorm(N)
t2 <- sqrt(rho)*t0 + sqrt(1-rho)*stats::rnorm(N)
dat1 <- sirt::sim.raschtype(t1 , b=seq(-1.5,1.5,length=I/2) ) 
dat2 <- sirt::sim.raschtype(t2 , b=seq(-1.5,1.5,length=I/2) ) 
dat <- as.matrix(cbind( dat1 , dat2) )
res <- sirt::unidim.test.csn( dat )
  ##  CSN Statistic = 0.06056 , p = 0.02
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Test for unidimensionality}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
