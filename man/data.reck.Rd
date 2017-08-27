%% File Name: data.reck.Rd
%% File Version: 0.20
%% File Last Change: 2017-01-18 18:08:37

\name{data.reck}

\alias{data.reck}
\alias{data.reck21}
\alias{data.reck61DAT1}
\alias{data.reck61DAT2}
\alias{data.reck73C1a}
\alias{data.reck73C1b}
\alias{data.reck75C2}
\alias{data.reck78ExA}
\alias{data.reck79ExB}

\docType{data}


\title{
Datasets from Reckase' Book \emph{Multidimensional Item Response Theory}
}
\description{
Some simulated datasets from Reckase (2009).
}
\usage{
data(data.reck21)
data(data.reck61DAT1)
data(data.reck61DAT2)
data(data.reck73C1a)
data(data.reck73C1b)
data(data.reck75C2)
data(data.reck78ExA)
data(data.reck79ExB)
}

\format{
\itemize{
\item The format of the \code{data.reck21} (Table 2.1, p. 45) is:

\code{List of 2} \cr
\code{ $ data: num [1:2500, 1:50] 0 0 0 1 1 0 0 0 1 0 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 2} \cr
\code{  .. ..$ : NULL} \cr
\code{  .. ..$ : chr [1:50] "I0001" "I0002" "I0003" "I0004" ...} \cr
\code{ $ pars:'data.frame':} \cr
\code{  ..$ a: num [1:50] 1.83 1.38 1.47 1.53 0.88 0.82 1.02 1.19 1.15 0.18 ...} \cr
\code{  ..$ b: num [1:50] 0.91 0.81 0.06 -0.8 0.24 0.99 1.23 -0.47 2.78 -3.85 ...} \cr
\code{  ..$ c: num [1:50] 0 0 0 0.25 0.21 0.29 0.26 0.19 0 0.21 ...} \cr

\item The format of the datasets \code{data.reck61DAT1} and
\code{data.reck61DAT2} (Table 6.1, p. 153) is

\code{List of 4} \cr
\code{ $ data : num [1:2500, 1:30] 1 0 0 1 1 0 0 1 1 0 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 2} \cr
\code{  .. ..$ : NULL} \cr
\code{  .. ..$ : chr [1:30] "A01" "A02" "A03" "A04" ...} \cr
\code{ $ pars :'data.frame':} \cr
\code{  ..$ a1: num [1:30] 0.747 0.46 0.861 1.014 0.552 ...} \cr
\code{  ..$ a2: num [1:30] 0.025 0.0097 0.0067 0.008 0.0204 0.0064 0.0861 ...} \cr
\code{  ..$ a3: num [1:30] 0.1428 0.0692 0.404 0.047 0.1482 ...} \cr
\code{  ..$ d : num [1:30] 0.183 -0.192 -0.466 -0.434 -0.443 ...} \cr
\code{ $ mu   : num [1:3] -0.4 -0.7 0.1} \cr
\code{ $ sigma: num [1:3, 1:3] 1.21 0.297 1.232 0.297 0.81 ...} \cr



The dataset \code{data.reck61DAT2} has correlated dimensions while
\code{data.reck61DAT1} has uncorrelated dimensions.

\item Datasets \code{data.reck73C1a} and \code{data.reck73C1b} use item parameters
from Table 7.3 (p. 188). The dataset \code{C1a} has uncorrelated dimensions,
while \code{C1b} has perfectly correlated dimensions. The items are sensitive to
3 dimensions. The format of the datasets is

\code{List of 4} \cr
\code{ $ data : num [1:2500, 1:30] 1 0 1 1 1 0 1 1 1 1 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 2} \cr
\code{  .. ..$ : NULL} \cr
\code{  .. ..$ : chr [1:30] "A01" "A02" "A03" "A04" ...} \cr
\code{ $ pars :'data.frame':  30 obs. of  4 variables:} \cr
\code{  ..$ a1: num [1:30] 0.747 0.46 0.861 1.014 0.552 ...} \cr
\code{  ..$ a2: num [1:30] 0.025 0.0097 0.0067 0.008 0.0204 0.0064 ...} \cr
\code{  ..$ a3: num [1:30] 0.1428 0.0692 0.404 0.047 0.1482 ...} \cr
\code{  ..$ d : num [1:30] 0.183 -0.192 -0.466 -0.434 -0.443 ...} \cr
\code{ $ mu   : num [1:3] 0 0 0} \cr
\code{ $ sigma: num [1:3, 1:3] 0.167 0.236 0.289 0.236 0.334 ...} \cr

\item The dataset \code{data.reck75C2} is simulated using item parameters
from Table 7.5 (p. 191). It contains items which are sensitive to only
one dimension but individuals which have abilities in three
uncorrelated dimensions. The format is

\code{List of 4} \cr
\code{ $ data : num [1:2500, 1:30] 0 0 1 1 1 0 0 1 1 1 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 2} \cr
\code{  .. ..$ : NULL} \cr
\code{  .. ..$ : chr [1:30] "A01" "A02" "A03" "A04" ...} \cr
\code{ $ pars :'data.frame':  30 obs. of  4 variables:} \cr
\code{  ..$ a1: num [1:30] 0.56 0.48 0.67 0.57 0.54 0.74 0.7 0.59 0.63 0.64 ...} \cr
\code{  ..$ a2: num [1:30] 0.62 0.53 0.63 0.69 0.58 0.69 0.75 0.63 0.64 0.64 ...} \cr
\code{  ..$ a3: num [1:30] 0.46 0.42 0.43 0.51 0.41 0.48 0.46 0.5 0.51 0.46 ...} \cr
\code{  ..$ d : num [1:30] 0.1 0.06 -0.38 0.46 0.14 0.31 0.06 -1.23 0.47 1.06 ...} \cr
\code{ $ mu   : num [1:3] 0 0 0} \cr
\code{ $ sigma: num [1:3, 1:3] 1 0 0 0 1 0 0 0 1} \cr

\item The dataset \code{data.reck78ExA} contains simulated item responses
from Table 7.8 (p. 204 ff.). There are three item clusters and
two ability dimensions. The format is

\code{List of 4} \cr
\code{ $ data : num [1:2500, 1:50] 0 1 1 0 1 0 0 0 0 0 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 2} \cr
\code{  .. ..$ : NULL} \cr
\code{  .. ..$ : chr [1:50] "A01" "A02" "A03" "A04" ...} \cr
\code{ $ pars :'data.frame':  50 obs. of  3 variables:} \cr
\code{  ..$ a1: num [1:50] 0.889 1.057 1.047 1.178 1.029 ...} \cr
\code{  ..$ a2: num [1:50] 0.1399 0.0432 0.016 0.0231 0.2347 ...} \cr
\code{  ..$ d : num [1:50] 0.2724 1.2335 -0.0918 -0.2372 0.8471 ...} \cr
\code{ $ mu   : num [1:2] 0 0} \cr
\code{ $ sigma: num [1:2, 1:2] 1 0 0 1} \cr

\item The dataset \code{data.reck79ExB} contains simulated item responses
from Table 7.9 (p. 207 ff.). There are three item clusters and
three ability dimensions. The format is

\code{List of 4} \cr
\code{ $ data : num [1:2500, 1:50] 1 1 0 1 0 0 0 1 1 0 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 2} \cr
\code{  .. ..$ : NULL} \cr
\code{  .. ..$ : chr [1:50] "A01" "A02" "A03" "A04" ...} \cr
\code{ $ pars :'data.frame':  50 obs. of  4 variables:} \cr
\code{  ..$ a1: num [1:50] 0.895 1.032 1.036 1.163 1.022 ...} \cr
\code{  ..$ a2: num [1:50] 0.052 0.132 0.144 0.13 0.165 ...} \cr
\code{  ..$ a3: num [1:50] 0.0722 0.1923 0.0482 0.1321 0.204 ...} \cr
\code{  ..$ d : num [1:50] 0.2724 1.2335 -0.0918 -0.2372 0.8471 ...} \cr
\code{ $ mu   : num [1:3] 0 0 0} \cr
\code{ $ sigma: num [1:3, 1:3] 1 0 0 0 1 0 0 0 1} \cr

}
  }
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Simulated datasets
}



\references{
Reckase, M. (2009). \emph{Multidimensional item response theory}. 
New York: Springer.
}

\examples{
\dontrun{	
#############################################################################
# EXAMPLE 1: data.reck21 dataset, Table 2.1, p. 45	
#############################################################################
data(data.reck21)

dat <- data.reck21$dat      # extract dataset

# items with zero guessing parameters
guess0 <- c( 1 , 2 , 3 , 9 ,11 ,27 ,30 ,35 ,45 ,49 ,50 )
I <- ncol(dat)

#***
# Model 1: 3PL estimation using rasch.mml2
est.c <- est.a <- 1:I
est.c[ guess0 ] <- 0
mod1 <- sirt::rasch.mml2( dat , est.a=est.a , est.c=est.c , mmliter= 300 )
summary(mod1)

#***
# Model 2: 3PL estimation using smirt
Q <- matrix(1,I,1)
mod2 <- sirt::smirt( dat , Qmatrix=Q , est.a= "2PL"  , est.c=est.c , increment.factor=1.01)
summary(mod2)

#***
# Model 3: estimation in mirt package
library(mirt)
itemtype <- rep("3PL" , I )
itemtype[ guess0 ] <- "2PL"
mod3 <- mirt::mirt(dat, 1, itemtype = itemtype , verbose=TRUE)
summary(mod3)

c3 <- unlist( coef(mod3) )[ 1:(4*I) ]
c3 <- matrix( c3 , I , 4 , byrow=TRUE )
# compare estimates of rasch.mml2, smirt and true parameters
round( cbind( mod1$item$c , mod2$item$c ,c3[,3]  ,data.reck21$pars$c ) , 2 )
round( cbind( mod1$item$a , mod2$item$a.Dim1 ,c3[,1], data.reck21$pars$a ) , 2 )
round( cbind( mod1$item$b , mod2$item$b.Dim1 / mod2$item$a.Dim1 , - c3[,2] / c3[,1] , 
            data.reck21$pars$b ) , 2 )
            
#############################################################################
# EXAMPLE 2: data.reck61 dataset, Table 6.1, p. 153	
#############################################################################            

data(data.reck61DAT1)
dat <- data.reck61DAT1$data

#***
# Model 1: Exploratory factor analysis

#-- Model 1a: tam.fa in TAM
library(TAM)
mod1a <- TAM::tam.fa( dat , irtmodel="efa" , nfactors=3 )
# varimax rotation
varimax(mod1a$B.stand)

# Model 1b: EFA in NOHARM (Promax rotation)
mod1b <- sirt::R2noharm( dat = dat , model.type="EFA" ,  dimensions = 3 , 
              writename = "reck61__3dim_efa", noharm.path = "c:/NOHARM" ,dec = ",")
summary(mod1b)

# Model 1c: EFA with noharm.sirt
mod1c <- sirt::noharm.sirt( dat=dat , dimensions=3  )
summary(mod1c)
plot(mod1c)

# Model 1d: EFA with 2 dimensions in noharm.sirt
mod1d <- sirt::noharm.sirt( dat=dat , dimensions=2  )
summary(mod1d)
plot(mod1d , efa.load.min=.2)   # plot loadings of at least .20

#***
# Model 2: Confirmatory factor analysis

#-- Model 2a: tam.fa in TAM
dims <- c( rep(1,10) , rep(3,10) , rep(2,10)  )
Qmatrix <- matrix( 0 , nrow=30 , ncol=3 )
Qmatrix[ cbind( 1:30 , dims) ] <- 1
mod2a <- TAM::tam.mml.2pl( dat ,Q=Qmatrix , 
            control=list( snodes=1000, QMC=TRUE , maxiter=200) )
summary(mod2a)

#-- Model 2b: smirt in sirt
mod2b <- sirt::smirt( dat ,Qmatrix =Qmatrix , est.a="2PL" , maxiter=20 , qmcnodes=1000 )
summary(mod2b)

#-- Model 2c: rasch.mml2 in sirt
mod2c <- sirt::rasch.mml2( dat ,Qmatrix =Qmatrix , est.a= 1:30 , 
                mmliter =200 , theta.k = seq(-5,5,len=11) )
summary(mod2c)

#-- Model 2d: mirt in mirt
cmodel <- mirt::mirt.model("
     F1 = 1-10
     F2 = 21-30
     F3 = 11-20
     COV = F1*F2, F1*F3 , F2*F3 " )
mod2d <- mirt::mirt(dat, cmodel , verbose=TRUE)
summary(mod2d)
coef(mod2d)

#-- Model 2e: CFA in NOHARM
# specify covariance pattern
P.pattern <- matrix( 1 , ncol=3 , nrow=3 )
P.init <- .4*P.pattern
diag(P.pattern) <- 0
diag(P.init) <- 1
# fix all entries in the loading matrix to 1
F.pattern <- matrix( 0 , nrow=30 , ncol=3 )
F.pattern[1:10,1] <- 1
F.pattern[21:30,2] <- 1
F.pattern[11:20,3] <- 1
F.init <- F.pattern
# estimate model
mod2e <- sirt::R2noharm( dat = dat , model.type="CFA" , P.pattern=P.pattern, 
            P.init=P.init , F.pattern=F.pattern, F.init=F.init ,
            writename = "reck61__3dim_cfa", noharm.path = "c:/NOHARM" ,dec = ",")
summary(mod2e)

#-- Model 2f: CFA with noharm.sirt
mod2f <- sirt::noharm.sirt( dat = dat , Fval=F.init , Fpatt = F.pattern ,
                 Pval=P.init , Ppatt = P.pattern )                      
summary(mod2f)

#############################################################################
# EXAMPLE 3: DETECT analysis data.reck78ExA and data.reck79ExB
############################################################################# 

data(data.reck78ExA)
data(data.reck79ExB)

#************************
# Example A
dat <- data.reck78ExA$data
#- estimate person score
score <- stats::qnorm( ( rowMeans( dat )+.5 )  / ( ncol(dat) + 1 ) )
#- extract item cluster
itemcluster <- substring( colnames(dat) , 1 , 1 )
#- confirmatory DETECT Item cluster
detectA <- sirt::conf.detect( data = dat , score = score , itemcluster = itemcluster )
  ##          unweighted weighted
  ##   DETECT      0.571    0.571
  ##   ASSI        0.523    0.523
  ##   RATIO       0.757    0.757

#- exploratory DETECT analysis
detect_explA <- sirt::expl.detect(data=dat, score, nclusters=10, N.est = nrow(dat)/2  )
  ##  Optimal Cluster Size is  5  (Maximum of DETECT Index)
  ##     N.Cluster N.items N.est N.val         size.cluster DETECT.est ASSI.est
  ##   1         2      50  1250  1250                31-19      0.531    0.404
  ##   2         3      50  1250  1250             10-19-21      0.554    0.407
  ##   3         4      50  1250  1250           10-19-14-7      0.630    0.509
  ##   4         5      50  1250  1250         10-19-3-7-11      0.653    0.546
  ##   5         6      50  1250  1250       10-12-7-3-7-11      0.593    0.458
  ##   6         7      50  1250  1250      10-12-7-3-7-9-2      0.604    0.474
  ##   7         8      50  1250  1250    10-12-7-3-3-9-4-2      0.608    0.481
  ##   8         9      50  1250  1250  10-12-7-3-3-5-4-2-4      0.617    0.494
  ##   9        10      50  1250  1250 10-5-7-7-3-3-5-4-2-4      0.592    0.460

# cluster membership
cluster_membership <- detect_explA$itemcluster$cluster3
# Cluster 1:
colnames(dat)[ cluster_membership == 1 ]
  ##   [1] "A01" "A02" "A03" "A04" "A05" "A06" "A07" "A08" "A09" "A10"
# Cluster 2:
colnames(dat)[ cluster_membership == 2 ]
  ##    [1] "B11" "B12" "B13" "B14" "B15" "B16" "B17" "B18" "B19" "B20" "B21" "B22"
  ##   [13] "B23" "B25" "B26" "B27" "B28" "B29" "B30"
# Cluster 3:
colnames(dat)[ cluster_membership == 3 ]
  ##    [1] "B24" "C31" "C32" "C33" "C34" "C35" "C36" "C37" "C38" "C39" "C40" "C41"
  ##   [13] "C42" "C43" "C44" "C45" "C46" "C47" "C48" "C49" "C50"

#************************
# Example B
dat <- data.reck79ExB$data
#- estimate person score
score <- stats::qnorm( ( rowMeans( dat )+.5 )  / ( ncol(dat) + 1 ) )
#- extract item cluster
itemcluster <- substring( colnames(dat) , 1 , 1 )
#- confirmatory DETECT Item cluster
detectB <- sirt::conf.detect( data = dat , score = score , itemcluster = itemcluster )
  ##          unweighted weighted
  ##   DETECT      0.715    0.715
  ##   ASSI        0.624    0.624
  ##   RATIO       0.855    0.855

#- exploratory DETECT analysis
detect_explB <- sirt::expl.detect(data=dat, score, nclusters=10, N.est = nrow(dat)/2  )
  ##   Optimal Cluster Size is  4  (Maximum of DETECT Index)
  ##   
  ##     N.Cluster N.items N.est N.val         size.cluster DETECT.est ASSI.est
  ##   1         2      50  1250  1250                30-20      0.665    0.546
  ##   2         3      50  1250  1250             10-20-20      0.686    0.585
  ##   3         4      50  1250  1250           10-20-8-12      0.728    0.644
  ##   4         5      50  1250  1250         10-6-14-8-12      0.654    0.553
  ##   5         6      50  1250  1250       10-6-14-3-12-5      0.659    0.561
  ##   6         7      50  1250  1250      10-6-14-3-7-5-5      0.664    0.576
  ##   7         8      50  1250  1250     10-6-7-7-3-7-5-5      0.616    0.518
  ##   8         9      50  1250  1250   10-6-7-7-3-5-5-5-2      0.612    0.512
  ##   9        10      50  1250  1250 10-6-7-7-3-5-3-5-2-2      0.613    0.512
}
}
\keyword{datasets}
