\name{wle.rasch.jackknife}
\alias{wle.rasch.jackknife}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standard Error Estimation of WLE by Jackknifing
}
\description{
This function calculates standard errors of WLEs (Warm, 1989) for 
stratified item designs and item designs with testlets for the 
Rasch model.
}
\usage{
wle.rasch.jackknife(dat, b, itemweights = 1 + 0 * b, pid = NULL, 
    testlet = NULL, stratum = NULL, size.itempop = NULL)
}
\arguments{
  \item{dat}{
An \eqn{N \times I} data frame of item responses
}
  \item{b}{
Vector of item difficulties
}
  \item{itemweights}{
Weights for items, i.e. fixed item discriminations
}
  \item{pid}{
Person identifier
}
  \item{testlet}{
A vector of length \eqn{I} which defines which item 
belongs to which testlet. If some items
does not belong to any testlet, then define separate
testlet labels for these single items.
}
  \item{stratum}{
Item stratum
}
  \item{size.itempop}{
Number of items in an item stratum of the finite item population.
}
}
\details{
The idea of Jackknife in item response models can be found in
Wainer and Wright (1980).
}
\value{
A list with following entries:
  \item{wle}{Data frame with some estimated statistics. The column
  \code{wle} is the WLE and \code{wle.jackse} its corresponding
  standard error estimated by jackknife.
  }
  \item{wle.rel}{WLE reliability (Adams, 2005)}
}
\references{
Adams, R. J. (2005). Reliability as a measurement design effect. 
\emph{Studies in Educational Evaluation}, \bold{31}, 162-172.

Gershunskaya, J., Jiang, J., & Lahiri, P. (2009). Resampling methods in surveys. 
In D. Pfeffermann and C.R. Rao (Eds.). \emph{Handbook of Statistics 29B; Sample
Surveys: Inference and Analysis} (pp. 121-151). Amsterdam: North Holland.

Wainer, H., & Wright, B. D. (1980). Robust estimation of ability in the Rasch 
model. \emph{Psychometrika}, \bold{45}, 373-391.

Warm, T. A. (1989). Weighted likelihood estimation of ability in item response theory.
\emph{Psychometrika}, \bold{54}, 427-450.
}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{wle.rasch}}
}
\examples{
#############################################################################
# EXAMPLE 1: Dataset Reading
#############################################################################
data(data.read)
dat <- data.read

# estimation of the Rasch model
res <- sirt::rasch.mml2( dat , parm.conv = .001)

# WLE estimation
wle1 <- sirt::wle.rasch(dat, b = res$item$thresh )

# simple jackknife WLE estimation
wle2 <- sirt::wle.rasch.jackknife(dat, b =res$item$thresh )
  ## WLE Reliability = 0.651 

# SE(WLE) for testlets A, B and C
wle3 <- sirt::wle.rasch.jackknife(dat, b =res$item$thresh , 
           testlet = substring( colnames(dat),1,1) )
  ## WLE Reliability = 0.572 
		
# SE(WLE) for item strata A,B, C
wle4 <- sirt::wle.rasch.jackknife(dat, b =res$item$thresh , 
             stratum = substring( colnames(dat),1,1) )
  ## WLE Reliability = 0.683 

# SE (WLE) for finite item strata
# A (10 items) , B (7 items) , C (4 items -> no sampling error)
# in every stratum 4 items were sampled
size.itempop <- c(10,7,4)
names(size.itempop) <- c("A","B","C")
wle5 <- sirt::wle.rasch.jackknife(dat, b =res$item$thresh , 
             stratum = substring( colnames(dat),1,1) ,
             size.itempop = size.itempop )
  ## Stratum  A (Mean) Correction Factor 0.6 
  ## Stratum  B (Mean) Correction Factor 0.42857 
  ## Stratum  C (Mean) Correction Factor 0 
  ## WLE Reliability = 0.876 
         
# compare different estimated standard errors
a2 <- stats::aggregate( wle2$wle$wle.jackse , list( wle2$wle$wle) , mean )
colnames(a2) <- c("wle" , "se.simple")
a2$se.testlet <- stats::aggregate( wle3$wle$wle.jackse , list( wle3$wle$wle) , mean )[,2]
a2$se.strata <- stats::aggregate( wle4$wle$wle.jackse , list( wle4$wle$wle) , mean )[,2]
a2$se.finitepop.strata <- stats::aggregate( wle5$wle$wle.jackse , 
	list( wle5$wle$wle) , mean )[,2]
round( a2 , 3 )
  ## > round( a2 , 3 )
  ##       wle se.simple se.testlet se.strata se.finitepop.strata
  ## 1  -5.085     0.440      0.649     0.331               0.138
  ## 2  -3.114     0.865      1.519     0.632               0.379
  ## 3  -2.585     0.790      0.849     0.751               0.495
  ## 4  -2.133     0.715      1.177     0.546               0.319
  ## 5  -1.721     0.597      0.767     0.527               0.317
  ## 6  -1.330     0.633      0.623     0.617               0.377
  ## 7  -0.942     0.631      0.643     0.604               0.365
  ## 8  -0.541     0.655      0.678     0.617               0.384
  ## 9  -0.104     0.671      0.646     0.659               0.434
  ## 10  0.406     0.771      0.706     0.751               0.461
  ## 11  1.080     1.118      0.893     1.076               0.630
  ## 12  2.332     0.400      0.631     0.272               0.195
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Weighted likelihood estimation (WLE)}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
