%% File Name: dmlavaan.Rd
%% File Version: 0.12

\name{dmlavaan}
\alias{dmlavaan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Comparing Regression Parameters of Different lavaan Models Fitted to the
Same Dataset
}
\description{
The function \code{dmlavaan} compares model parameters from different \pkg{lavaan}
models fitted to the same dataset. This leads to dependent coefficients.
Statistical inference is either conducted by M-estimation (i.e., robust
sandwich method; \code{method="bootstrap"}) or bootstrap (\code{method="bootstrap"}).
See Mize et al. (2019) or Weesie (1999) for more details.
}
\usage{
dmlavaan(fun1, args1, fun2, args2, method="sandwich", R=50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fun1}{
lavaan function of the first model (e.g., \code{"lavaan"}, \code{"cfa"}, or
\code{"sem"})
}
  \item{args1}{
arguments for lavaan function in the first model
}
  \item{fun2}{
lavaan function of the second model (e.g., \code{"lavaan"}, \code{"cfa"}, or
\code{"sem"})
}
  \item{args2}{
arguments for lavaan function in the second model
}
  \item{method}{
estimation method for standard errors
}
  \item{R}{
Number of bootstrap samples
}
}

\details{
In bootstrap estimation, a normal approximation is applied in the
computation of confidence intervals. Hence, \code{R} could be chosen
relatively small. \cr \cr
TO DO (not yet implemented): \cr
\tabular{ll}{
1) \tab inclusion of sampling weights \cr
2) \tab cluster robust standard errors in hierarchical sampling \cr
3) \tab stratification \cr
}
}

\value{
A list with following entries
\item{coef}{Model parameters of both models}
\item{vcov}{Covariance matrix of model parameters of both models}
\item{partable}{Parameter table containing all univariate model parameters}
\item{\ldots}{More entries}
}

\references{
Mize, T.D., Doan, L., & Long, J.S. (2019).
A general framework for comparing predictions and marginal effects across models.
\emph{Sociological Methodology, 49}(1), 152-189.
\doi{10.1177/0081175019852763}

Weesie, J. (1999) Seemingly unrelated estimation and the cluster-adjusted sandwich
estimator. \emph{Stata Technical Bulletin, 9}, 231-248.
}


%\author{
%%  ~~who you are~~
%}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}


\examples{
\dontrun{
############################################################################
# EXAMPLE 1: Confirmatory factor analysis with and without fourth item
#############################################################################

#**** simulate data
N <- 200  # number of persons
I <- 4    # number of items

# loadings and error correlations
lam <- seq(.7,.4, len=I)
PSI <- diag( 1-lam^2 )

# define some model misspecification
sd_error <- .1
S1 <- matrix( c( -1.84, 0.39,-0.68, 0.13,
  0.39,-1.31,-0.07,-0.27,
 -0.68,-0.07, 0.90, 1.91,
  0.13,-0.27, 1.91,-0.56 ), nrow=4, ncol=4, byrow=TRUE)
S1 <- ( S1 - mean(S1) ) / sd(S1) * sd_error

Sigma <- lam \%*\% t(lam) + PSI + S1
dat <- MASS::mvrnorm(n=N, mu=rep(0,I), Sigma=Sigma)
colnames(dat) <- paste0("X",1:4)
dat <- as.data.frame(dat)
rownames(Sigma) <- colnames(Sigma) <- colnames(dat)


#*** define two lavaan models
lavmodel1 <- "F=~ X1 + X2 + X3 + X4"
lavmodel2 <- "F=~ X1 + X2 + X3"

#*** define lavaan estimation arguments and functions
fun2 <- fun1 <- "cfa"
args1 <- list( model=lavmodel1, data=dat, std.lv=TRUE, estimator="MLR")
args2 <- args1
args2$model <- lavmodel2

#* run model comparison
res1 <- sirt::dmlavaan( fun1=fun1, args1=args1, fun2=fun2, args2=args2)

# inspect results
sirt:::print_digits(res1$partable, digits=3)
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
