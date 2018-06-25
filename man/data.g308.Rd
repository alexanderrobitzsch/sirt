%% File Name: data.g308.Rd
%% File Version: 0.15

\name{data.g308}
\alias{data.g308}

\docType{data}
\title{
C-Test Datasets
}

\description{
Some datasets of C-tests are provided. The dataset \code{data.g308}
was used in Schroeders, Robitzsch and Schipolowski (2014).
}

\usage{
data(data.g308)
}

\format{
\itemize{
\item The dataset \code{data.g308} is a C-test containing 20 items and is
used in Schroeders, Robitzsch and Schipolowski (2014) and is of the
following format \cr

\code{'data.frame':   747 obs. of  21 variables:} \cr
\code{ $ id    : int  1 2 3 4 5 6 7 8 9 10 ...} \cr
\code{ $ G30801: int  1 1 1 1 1 0 0 1 1 1 ...} \cr
\code{ $ G30802: int  1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{ $ G30803: int  1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{ $ G30804: int  1 1 1 1 1 0 1 1 1 1 ...} \cr
\code{[...]} \cr
\code{ $ G30817: int  0 0 0 0 1 0 1 0 1 0 ...} \cr
\code{ $ G30818: int  0 0 1 0 0 0 0 1 1 0 ...} \cr
\code{ $ G30819: int  1 1 1 1 0 0 1 1 1 0 ...} \cr
\code{ $ G30820: int  1 1 1 1 0 0 0 1 1 0 ...} \cr
}
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}

\references{
Schroeders, U., Robitzsch, A., & Schipolowski, S. (2014).
A comparison of different psychometric approaches to modeling testlet
structures: An example with C-tests.
\emph{Journal of Educational Measurement}, \bold{51(4)}, 400-418.
}

\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Dataset G308 from Schroeders et al. (2014)
#############################################################################

data(data.g308)
dat <- data.g308

library(TAM)
library(sirt)
library(combinat)

# define testlets
testlet <- c(1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 4, 5, 5, 6, 6, 6)

#****************************************
#*** Model 1: Rasch model
mod1 <- TAM::tam.mml(resp=dat[,-1], pid=dat[,1],
       control=list(maxiter=300, snodes=1500))
summary(mod1)

#****************************************
#*** Model 2: Rasch testlet model

# testlets are dimensions, assign items to Q-matrix
TT <- length(unique(testlet))
Q <- matrix(0, nrow=ncol(dat)-1, ncol=TT + 1)
Q[,1] <- 1 # First dimension constitutes g-factor
for (tt in 1:TT){Q[testlet==tt, tt+1] <- 1}

# In a testlet model, all dimensions are uncorrelated among
# each other, that is, all pairwise correlations are set to 0,
# which can be accomplished with the "variance.fixed" command
variance.fixed <- cbind(t( combinat::combn(TT+1,2)), 0)
mod2 <- TAM::tam.mml(resp=dat[,-1], pid=dat[,1], Q=Q,
            variance.fixed=variance.fixed,
            control=list(snodes=1500, maxiter=300))
summary(mod2)

#****************************************
#*** Model 3: Partial credit model

scores <- list()
testlet.names <- NULL
dat.pcm <- NULL
for (tt in 1:max(testlet) ){
   scores[[tt]] <- rowSums (dat[,-1][, testlet==tt, drop=FALSE])
   dat.pcm <- c(dat.pcm, list(c(scores[[tt]])))
   testlet.names <- append(testlet.names, paste0("testlet",tt) )
   }
dat.pcm <- as.data.frame(dat.pcm)
colnames(dat.pcm) <- testlet.names
mod3 <- TAM::tam.mml(resp=dat.pcm, control=list(snodes=1500, maxiter=300) )
summary(mod3)

#****************************************
#*** Model 4: Copula model

mod4 <- sirt::rasch.copula2 (dat=dat[,-1], itemcluster=testlet)
summary(mod4)
}
}



\keyword{datasets}
