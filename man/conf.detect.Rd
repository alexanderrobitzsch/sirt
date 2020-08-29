%% File Name: conf.detect.Rd
%% File Version: 1.479

\name{conf.detect}
\alias{conf.detect}
\alias{summary.conf.detect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confirmatory DETECT and polyDETECT Analysis
}
\description{
This function computes the DETECT statistics for dichotomous item responses
and the polyDETECT statistic for polytomous item responses
under a confirmatory specification of item clusters
(Stout, Habing, Douglas & Kim, 1996; Zhang & Stout,
1999a, 1999b; Zhang, 2007; Bonifay, Reise, Scheines, & Meijer, 2015).

Item responses in a multi-matrix design are allowed (Zhang, 2013).

An exploratory DETECT analysis can be conducted using the
\code{\link{expl.detect}} function.
}
\usage{
conf.detect(data, score, itemcluster, bwscale=1.1, progress=TRUE,
        thetagrid=seq(-3, 3, len=200), smooth=TRUE, use_sum_score=FALSE)

\method{summary}{conf.detect}(object, digits=3, file=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
An \eqn{N \times I} data frame of dichotomous or polytomous responses.
Missing responses are allowed.
}
  \item{score}{
An ability estimate, e.g. the WLE, sum score or mean score
}
  \item{itemcluster}{
Item cluster for each item. The order of entries must correspond
to the columns in \code{data}.
}
  \item{bwscale}{
Bandwidth factor for calculation of conditional covariance
(see \code{\link{ccov.np}})
}
  \item{progress}{
Display progress?
}
\item{smooth}{Logical indicating whether smoothing should be
applied for conditional covariance estimation}
  \item{thetagrid}{
A vector which contains theta values where conditional
covariances are evaluated.
}
\item{use_sum_score}{Logical indicating whether sum score should be used.
With this option, the bias corrected conditional covariance of Zhang and
Stout (1999) is used.}
\item{object}{Object of class \code{conf.detect}}
\item{digits}{Number of digits for rounding in \code{summary}}
\item{file}{Optional file name to be sunk for \code{summary}}
\item{\dots}{Further arguments to be passed}
}
\details{
The result of DETECT are the indices \code{DETECT}, \code{ASSI}
and \code{RATIO} (see Zhang 2007 for details) calculated
for the options \code{unweighted} and \code{weighted}.
The option \code{unweighted} means that all conditional covariances of
item pairs are equally weighted, \code{weighted} means that
these covariances are weighted by the sample size of item pairs.
In case of multi matrix item designs, both types of indices can
differ.

The classification scheme of these indices are as follows
(Jang & Roussos, 2007; Zhang, 2007):

\tabular{ll}{
 Strong multidimensionality \tab DETECT > 1.00 \cr
Moderate multidimensionality \tab .40 < DETECT < 1.00 \cr
Weak multidimensionality \tab .20 < DETECT < .40 \cr
Essential unidimensionality \tab DETECT < .20 \cr
    }
\tabular{lll}{
 Maximum value under simple structure \tab ASSI=1 \tab RATIO=1 \cr
Essential deviation from unidimensionality \tab ASSI > .25 \tab
RATIO > .36 \cr
Essential unidimensionality \tab  ASSI < .25 \tab
RATIO < .36 \cr
    }


Note that the expected value of a conditional covariance for an item pair
is negative when a unidimensional model holds. In consequence,
the DETECT index can become negative for unidimensional data
(see Example 3). This can be also seen in the statistic
\code{MCOV100} in the value \code{detect}.
}
\value{
A list with following entries:
\item{detect}{Data frame with statistics DETECT, ASSI, RATIO, MADCOV100
and MCOV100}
\item{ccovtable}{Individual contributions to conditional covariance}
\item{ccov.matrix}{Evaluated conditional covariance}
}


\references{
Bonifay, W. E., Reise, S. P., Scheines, R., & Meijer, R. R. (2015).
When are multidimensional data unidimensional enough for structural
equation modeling? An evaluation of the DETECT multidimensionality index.
\emph{Structural Equation Modeling, 22}(4), 504-516.
\doi{10.1080/10705511.2014.938596}

Jang, E. E., & Roussos, L. (2007). An investigation into the dimensionality
of TOEFL using conditional covariance-based nonparametric approach.
\emph{Journal of Educational Measurement, 44}(1), 1-21.
\doi{10.1111/j.1745-3984.2007.00024.x}

Stout, W., Habing, B., Douglas, J., & Kim, H. R. (1996).
Conditional covariance-based nonparametric multidimensionality assessment.
\emph{Applied Psychological Measurement, 20}(4), 331-354.
\doi{10.1177/014662169602000403}

Zhang, J. (2007). Conditional covariance theory and DETECT for
polytomous items. \emph{Psychometrika, 72}(1), 69-91.
\doi{10.1007/s11336-004-1257-7}

Zhang, J. (2013). A procedure for dimensionality analyses of
response data from various test designs. \emph{Psychometrika, 78}(1), 37-58.
\doi{10.1007/s11336-012-9287-z}

Zhang, J., & Stout, W. (1999a). Conditional covariance structure
of generalized compensatory multidimensional items.
\emph{Psychometrika, 64}(2), 129-152.
\doi{10.1007/BF02294532}

Zhang, J., & Stout, W. (1999b). The theoretical DETECT index of
dimensionality and its application to approximate simple structure.
\emph{Psychometrika, 64}(2), 213-249.
\doi{10.1007/BF02294536}
}


%\author{
%Alexander Robitzsch
%}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For a download of the free \emph{DIM-Pack} software (DIMTEST, DETECT) see
\emph{https://psychometrics.onlinehelp.measuredprogress.org/tools/dim/}.

See \code{\link{expl.detect}} for exploratory DETECT analysis.
}


\examples{
#############################################################################
# EXAMPLE 1: TIMSS mathematics data set (dichotomous data)
#############################################################################
data(data.timss)

# extract data
dat <- data.timss$data
dat <- dat[, substring( colnames(dat),1,1)=="M" ]
# extract item informations
iteminfo <- data.timss$item
# estimate Rasch model
mod1 <- sirt::rasch.mml2( dat )
# estimate WLEs
wle1 <- sirt::wle.rasch( dat, b=mod1$item$b )$theta

# DETECT for content domains
detect1 <- sirt::conf.detect( data=dat, score=wle1,
                    itemcluster=iteminfo$Content.Domain )
  ##          unweighted weighted
  ##   DETECT      0.316    0.316
  ##   ASSI        0.273    0.273
  ##   RATIO       0.355    0.355

\dontrun{
# DETECT cognitive domains
detect2 <- sirt::conf.detect( data=dat, score=wle1,
                    itemcluster=iteminfo$Cognitive.Domain )
  ##          unweighted weighted
  ##   DETECT      0.251    0.251
  ##   ASSI        0.227    0.227
  ##   RATIO       0.282    0.282

# DETECT for item format
detect3 <- sirt::conf.detect( data=dat, score=wle1,
                    itemcluster=iteminfo$Format )
  ##          unweighted weighted
  ##   DETECT      0.056    0.056
  ##   ASSI        0.060    0.060
  ##   RATIO       0.062    0.062

# DETECT for item blocks
detect4 <- sirt::conf.detect( data=dat, score=wle1,
                    itemcluster=iteminfo$Block )
  ##          unweighted weighted
  ##   DETECT      0.301    0.301
  ##   ASSI        0.193    0.193
  ##   RATIO       0.339    0.339 }

# Exploratory DETECT: Application of a cluster analysis employing the Ward method
detect5 <- sirt::expl.detect( data=dat, score=wle1,
                nclusters=10, N.est=nrow(dat)  )
# Plot cluster solution
pl <- graphics::plot( detect5$clusterfit, main="Cluster solution" )
stats::rect.hclust(detect5$clusterfit, k=4, border="red")

\dontrun{
#############################################################################
# EXAMPLE 2: Big 5 data set (polytomous data)
#############################################################################

# attach Big5 Dataset
data(data.big5)

# select 6 items of each dimension
dat <- data.big5
dat <- dat[, 1:30]

# estimate person score by simply using a transformed sum score
score <- stats::qnorm( ( rowMeans( dat )+.5 )  / ( 30 + 1 ) )

# extract item cluster (Big 5 dimensions)
itemcluster <- substring( colnames(dat), 1, 1 )

# DETECT Item cluster
detect1 <- sirt::conf.detect( data=dat, score=score, itemcluster=itemcluster )
  ##        unweighted weighted
  ## DETECT      1.256    1.256
  ## ASSI        0.384    0.384
  ## RATIO       0.597    0.597

# Exploratory DETECT
detect5 <- sirt::expl.detect( data=dat, score=score,
                     nclusters=9, N.est=nrow(dat)  )
  ## DETECT (unweighted)
  ## Optimal Cluster Size is  6  (Maximum of DETECT Index)
  ##   N.Cluster N.items N.est N.val      size.cluster DETECT.est ASSI.est RATIO.est
  ## 1         2      30   500     0              6-24      1.073    0.246     0.510
  ## 2         3      30   500     0           6-10-14      1.578    0.457     0.750
  ## 3         4      30   500     0         6-10-11-3      1.532    0.444     0.729
  ## 4         5      30   500     0        6-8-11-2-3      1.591    0.462     0.757
  ## 5         6      30   500     0       6-8-6-2-5-3      1.610    0.499     0.766
  ## 6         7      30   500     0     6-3-6-2-5-5-3      1.557    0.476     0.740
  ## 7         8      30   500     0   6-3-3-2-3-5-5-3      1.540    0.462     0.732
  ## 8         9      30   500     0 6-3-3-2-3-5-3-3-2      1.522    0.444     0.724

# Plot Cluster solution
pl <- graphics::plot( detect5$clusterfit, main="Cluster solution" )
stats::rect.hclust(detect5$clusterfit, k=6, border="red")

#############################################################################
# EXAMPLE 3: DETECT index for unidimensional data
#############################################################################

set.seed(976)
N <- 1000
I <- 20
b <- sample( seq( -2, 2, len=I) )
dat <- sirt::sim.raschtype( stats::rnorm(N), b=b )

# estimate Rasch model and corresponding WLEs
mod1 <- TAM::tam.mml( dat )
wmod1 <- TAM::tam.wle(mod1)$theta

# define item cluster
itemcluster <- c( rep(1,5), rep(2,I-5) )

# compute DETECT statistic
detect1 <- sirt::conf.detect( data=dat, score=wmod1, itemcluster=itemcluster)
  ##            unweighted weighted
  ##  DETECT        -0.184   -0.184
  ##  ASSI          -0.147   -0.147
  ##  RATIO         -0.226   -0.226
  ##  MADCOV100      0.816    0.816
  ##  MCOV100       -0.786   -0.786
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{DETECT}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
