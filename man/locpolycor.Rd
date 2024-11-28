%% File Name: locpolycor.Rd
%% File Version: 0.07

\name{locpolycor}
\alias{locpolycor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Local Modeling of Thresholds and Polychoric Correlations
}
\description{
Estimates thresholds and polychoric correlations as a function of a
continuous moderator variable \eqn{x}.
}
\usage{
locpolycor(y, data.mod, moderator.grid, h=1.1, model_thresh, model_polycor,
     kernel="gaussian", eps=1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Matrix with columns referrring to ordinal items
}
  \item{data.mod}{
Values of the moderator variable \eqn{x}
}
  \item{moderator.grid}{
Grid of \eqn{x} values to be used for local estimation of thresholds
and polychoric correlations
}
  \item{h}{
Bandwidth factor
}
  \item{model_thresh}{
Model for thresholds: can be \code{'const'} (constant function) or \code{'lin'}
(linear function)
}
  \item{model_polycor}{
Model for polychoric correlations: can be \code{'const'} (constant function) or
\code{'lin'} (linear function)
}
  \item{kernel}{
Used kernel function, see \code{\link{lsem_local_weights}}
}
  \item{eps}{
Parameter added in likelihood optimization
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
A list with entries
\item{thresh_list}{Threshold parameters}
\item{thresh_stat}{Estimated thresholds}
\item{polycor_stat}{Estimated polychoric correlations}
\item{\ldots}{...}
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}


\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Two items, moderator on (0,1)
#############################################################################


#*** simulate data

# functions for thresholds
th1_fun <- function(x){ -0.3*(x-2)^2 + .2 }
th2_fun <- function(x){ 0.4*(x+1)^2 - 0.6 }
zh1_fun <- function(x){ 0.3*(x-1) }
# function polychoric correlation
cor_x12 <- function(x){ 0.2+0.1*(x-0.5)+0.09*(x-0.5)^2 }

# simulate moderator
x <- stats::runif(N)

# simulate data
yast <- matrix( NA, nrow=N, ncol=2)
for (nn in 1:N){
    rho12 <- cor_x12(x[nn])
    Sigma <- matrix(0, 2,2)
    Sigma[1,2] <- rho12
    Sigma <- Sigma + t(Sigma)
    diag(Sigma) <- 1
    yast_nn <- MASS::mvrnorm( 1, mu=rep(0,2), Sigma=Sigma )
    yast[nn,] <- yast_nn
}
y <- 0*yast
th1_x <- th1_fun(x)
th2_x <- th2_fun(x)
zh1_x <- zh1_fun(x)
y[,1] <- 1*( yast[,1] > th1_x ) + 1*( yast[,1] > th2_x )
y[,2] <- 1*( yast[,2] > zh1_x )
colnames(y) <- paste0("I",1:ncol(y))

dat <- data.frame(x=x, y)

#-- local modeling
res <- sirt::locpolycor(y, data.mod=x, moderator.grid=c(0, .25, .5, .75, 1 ), h=2,
                 model_thresh="lin", model_polycor="lin")
str(res)
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
