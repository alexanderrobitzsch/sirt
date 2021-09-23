%% File Name: data.trees.Rd
%% File Version: 0.14

\name{data.trees}
\alias{data.trees}
\docType{data}

\title{
Dataset Used in Stoyan, Pommerening and Wuensche (2018)
}
\description{
Dataset used in Stoyan, Pommerening and Wuensche (2018; see also
Pommerening et al., 2018). In the dataset, 15 forest managers classify 387 trees
either as trees to be maintained or as trees to be removed. They assign tree marks,
either 0 or 1, where mark 1 means remove.
}

\usage{
data(data.trees)
}

\format{
The dataset has the following structure. \cr

\code{'data.frame':   387 obs. of  16 variables:} \cr
\code{ $ Number: int  142 184 9 300 374 42 382 108 125 201 ...} \cr
\code{ $ FM1   : int  1 1 1 1 1 1 1 1 1 0 ...} \cr
\code{ $ FM2   : int  1 1 1 0 1 1 1 1 1 1 ...} \cr
\code{ $ FM3   : int  1 0 1 1 1 1 1 1 1 1 ...} \cr
\code{ $ FM4   : int  1 1 1 1 1 1 0 1 1 1 ...} \cr
\code{ $ FM5   : int  1 1 1 1 1 1 0 0 0 1 ...} \cr
\code{ $ FM6   : int  1 1 1 1 0 1 1 1 1 0 ...} \cr
\code{ $ FM7   : int  1 0 1 1 0 0 1 0 1 1 ...} \cr
\code{ $ FM8   : int  1 1 1 1 1 0 0 1 0 1 ...} \cr
\code{ $ FM9   : int  1 1 0 1 1 1 1 0 1 1 ...} \cr
\code{ $ FM10  : int  0 1 1 0 1 1 1 1 0 0 ...} \cr
\code{ $ FM11  : int  1 1 1 1 0 1 1 0 1 0 ...} \cr
\code{ $ FM12  : int  1 1 1 1 1 1 0 1 0 0 ...} \cr
\code{ $ FM13  : int  0 1 0 0 1 1 1 1 1 1 ...} \cr
\code{ $ FM14  : int  1 1 1 1 1 0 1 1 1 1 ...} \cr
\code{ $ FM15  : int  1 1 0 1 1 0 1 0 0 1 ...} \cr

}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
\url{https://www.pommerening.org/wiki/images/d/dc/CoedyBreninSortedforPublication.txt}
}

\references{
Pommerening, A., Ramos, C. P., Kedziora, W., Haufe, J., & Stoyan, D. (2018).
Rating experiments in forestry: How much agreement is there in tree marking?
\emph{PloS ONE, 13}(3), e0194747. \doi{10.1371/journal.pone.0194747}

Stoyan, D., Pommerening, A., & Wuensche, A. (2018).
Rater classification by means of set-theoretic methods applied to forestry data.
\emph{Journal of Environmental Statistics, 8}(2), 1-17.
\url{http://www.jenvstat.org/v08/i02}
}


\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Latent class models, latent trait models, mixed membership models
#############################################################################

data(data.trees, package="sirt")
dat <- data.trees[,-1]
I <- ncol(dat)

#** latent class models with 2, 3, and 4 classes
problevels <- seq( 0, 1, len=2 )
mod02 <- sirt::gom.em(dat, K=2, problevels, model="GOM")
mod03 <- sirt::gom.em(dat, K=3, problevels, model="GOM")
mod04 <- sirt::gom.em(dat, K=4, problevels, model="GOM")

#** grade of membership models
mod11 <- sirt::gom.em(dat, K=2, theta0.k=10*seq(-1,1,len=11), model="GOMnormal")
problevels <- seq( 0, 1, len=3 )
mod12 <- sirt::gom.em(dat, K=2, problevels, model="GOM")
mod13 <- sirt::gom.em(dat, K=3, problevels, model="GOM")
mod14 <- sirt::gom.em(dat, K=4, problevels, model="GOM")
problevels <- seq( 0, 1, len=4 )
mod22 <- sirt::gom.em(dat, K=2, problevels, model="GOM")
mod23 <- sirt::gom.em(dat, K=3, problevels, model="GOM")
mod24 <- sirt::gom.em(dat, K=4, problevels, model="GOM")

#** latent trait models
#- 1PL
mod31 <- sirt::rasch.mml2(dat)
#- 2PL
mod32 <- sirt::rasch.mml2(dat, est.a=1:I)

#- model comparison
IRT.compareModels(mod02, mod03, mod04, mod11, mod12, mod13, mod14,
                     mod22, mod23, mod24, mod31, mod32)

#-- inspect model results
summary(mod12)
round( cbind( mod12$theta.k, mod12$pi.k ),3)

summary(mod13)
round(cbind( mod13$theta.k, mod13$pi.k ),3)
}
}



%% \keyword{datasets}
