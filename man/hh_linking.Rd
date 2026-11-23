%% File Name: hh_linking.Rd
%% File Version: 0.079

\name{hh_linking}
\alias{hh_linking}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Haberman-Halpin Linking for Multiple Groups
}
\description{
Implements Haberman-Halpin linking which is defined as pairwise Haberman linking
with the extensions for using precision weights and robust weights for outlying DIF
effects. The arguments can be differently chosen for item discriminations and
item intercepts.
}
\usage{
hh_linking(mirt_mg, prec_weights="none", rob_weights="none", eps=0.001,
    c=NULL, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mirt_mg}{
Object of class \code{\link{mirt_mg}}
}
  \item{prec_weights}{
Precision weights with options \code{"none"} or \code{"raw"} (precision weights).
Can also be a vector of length two if different specifications for item discriminations
and intercepts are requested
}
  \item{rob_weights}{
Type of robust weights. Can be \code{"bisquare1"} (significance-based approach
of Halpin), \code{"bisquare2"} (data-driven tuning constant in bisquare regression),
\code{"L0"} (\eqn{L_0} loss function) or \code{"Lp"}) (\eqn{L_p} loss function
with concrete values of \code{p})
}
  \item{eps}{
Tuning parameter for \code{L0} pr \code{Lp} loss function
}
  \item{c}{
Tuning value (constant factor) for \code{bisquare1} or \code{bisquare2}. If no values
are provided the defaults 1.96 and 3.0 are used, respectively.
}
\item{verbose}{Verbose argument}
}

\details{
The function implements pairwise Haberman linking with log-transformed item
discriminations and uses item intercepts instead of item difficulties.
}


\value{
List with following entries
\item{delta}{Estimated distribution parameters}
\item{res_sigma}{Informations about model fit for \code{sigma}}
\item{res_mu}{Informations about model fit for \code{mu}}
}

%\references{
%% ~put references to the literature/web site here ~
%}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mirt_mg}}, \code{\link{linking.haberman.lq}}
}


\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Simulated example with three groups
#############################################################################

set.seed(20260519)
I <- 10  # number of items
# simulate item parameters
pars1 <- data.frame( a=stats::runif(I, 0.5, 1.5), b=stats::runif(I, -2,2) )

itemprop <- .7  # proportion of administered items per group
N <- 500    # sample size per group

# distribution parameters
mu <- c(0,.3,-0.5)
sig <- c(1,1.2, 0.8)
G <- length(mu)

#--- simulate item responses
datasets <- list()
data <- NULL
for (gg in 1L:G){
  b1 <- pars1$b
  dat1 <- sirt::sim.raschtype( stats::rnorm(N, mean=mu[gg], sd=sig[gg]),
                                 b=b1, fixed.a=pars1$a )
  if (itemprop > 0){
    sel <- sort( sample( colnames(dat1), size=itemprop*I ) )
    dat1 <- dat1[,sel]
  }
  data1 <- data.frame( group=gg, dat1 )
  data <- plyr::rbind.fill( data, data1)
}

data <- data[, c("group", sort(colnames(data)[-1]) ) ]

#- groupwise 2PL fit
rmod <- sirt::mirt_mg(resp=data[,-1], group=data$group)

#- Haberman linking with item intercepts
mod0 <- sirt::hh_linking(rmod, prec_weights="none", rob_weights="none")
mod0$delta

#- Halpin approach using precision weights and bidquare loss function
mod1 <- sirt::hh_linking(rmod, prec_weights="raw", rob_weights="bisquare1")
mod1$delta
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
