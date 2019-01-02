%% File Name: data.long.Rd
%% File Version: 0.29

\name{data.long}
\alias{data.long}
\docType{data}
\title{
Longitudinal Dataset
}
\description{
This dataset contains 200 observations on
12 items. 6 items (\code{I1T1}, ...,\code{I6T1})
were administered at measurement occasion T1
and 6 items at T2 (\code{I3T2}, ..., \code{I8T2}). There were 4 anchor items
which were presented at both time points.
The first column in the dataset contains the student identifier.
}
\usage{
data(data.long)
}

\format{
The format of the dataset is

\code{'data.frame':   200 obs. of  13 variables:} \cr
\code{ $ idstud: int  1001 1002 1003 1004 1005 1006 1007 1008 1009 1010 ...} \cr
\code{ $ I1T1  : int  1 1 1 1 1 1 1 0 1 1 ...} \cr
\code{ $ I2T1  : int  0 0 1 1 1 1 0 1 1 1 ...} \cr
\code{ $ I3T1  : int  1 0 1 1 0 1 0 0 0 0 ...} \cr
\code{ $ I4T1  : int  1 0 0 1 0 0 0 0 1 1 ...} \cr
\code{ $ I5T1  : int  1 0 0 1 0 0 0 0 1 0 ...} \cr
\code{ $ I6T1  : int  1 0 0 0 0 0 0 0 0 0 ...} \cr
\code{ $ I3T2  : int  1 1 0 0 1 1 1 1 0 1 ...} \cr
\code{ $ I4T2  : int  1 1 0 0 1 1 0 0 0 1 ...} \cr
\code{ $ I5T2  : int  1 0 1 1 1 1 1 0 1 1 ...} \cr
\code{ $ I6T2  : int  1 1 0 0 0 0 0 0 0 1 ...} \cr
\code{ $ I7T2  : int  1 0 0 0 0 0 0 0 0 1 ...} \cr
\code{ $ I8T2  : int  0 0 0 0 1 0 0 0 0 0 ...} \cr
}





\examples{
\dontrun{
data(data.long)
dat <- data.long
dat <- dat[,-1]
I <- ncol(dat)

#*************************************************
# Model 1: 2-dimensional Rasch model
#*************************************************
# define Q-matrix
Q <- matrix(0,I,2)
Q[1:6,1] <- 1
Q[7:12,2] <- 1
rownames(Q) <- colnames(dat)
colnames(Q) <- c("T1","T2")

# vector with same items
itemnr <- as.numeric( substring( colnames(dat),2,2) )
# fix mean at T2 to zero
mu.fixed <- cbind( 2,0 )

#--- M1a: rasch.mml2 (in sirt)
mod1a <- sirt::rasch.mml2(dat, Q=Q, est.b=itemnr, mu.fixed=mu.fixed)
summary(mod1a)

#--- M1b: smirt (in sirt)
mod1b <- sirt::smirt(dat, Qmatrix=Q, irtmodel="comp", est.b=itemnr,
                  mu.fixed=mu.fixed )

#--- M1c: tam.mml (in TAM)

# assume equal item difficulty of I3T1 and I3T2, I4T1 and I4T2, ...
# create draft design matrix and modify it
A <- TAM::designMatrices(resp=dat)$A
dimnames(A)[[1]] <- colnames(dat)
  ##   > str(A)
  ##    num [1:12, 1:2, 1:12] 0 0 0 0 0 0 0 0 0 0 ...
  ##    - attr(*, "dimnames")=List of 3
  ##     ..$ : chr [1:12] "Item01" "Item02" "Item03" "Item04" ...
  ##     ..$ : chr [1:2] "Category0" "Category1"
  ##     ..$ : chr [1:12] "I1T1" "I2T1" "I3T1" "I4T1" ...
A1 <- A[,, c(1:6, 11:12 ) ]
A1[7,2,3] <- -1     # difficulty(I3T1)=difficulty(I3T2)
A1[8,2,4] <- -1     # I4T1=I4T2
A1[9,2,5] <- A1[10,2,6] <- -1
dimnames(A1)[[3]] <- substring( dimnames(A1)[[3]],1,2)
  ##   > A1[,2,]
  ##        I1 I2 I3 I4 I5 I6 I7 I8
  ##   I1T1 -1  0  0  0  0  0  0  0
  ##   I2T1  0 -1  0  0  0  0  0  0
  ##   I3T1  0  0 -1  0  0  0  0  0
  ##   I4T1  0  0  0 -1  0  0  0  0
  ##   I5T1  0  0  0  0 -1  0  0  0
  ##   I6T1  0  0  0  0  0 -1  0  0
  ##   I3T2  0  0 -1  0  0  0  0  0
  ##   I4T2  0  0  0 -1  0  0  0  0
  ##   I5T2  0  0  0  0 -1  0  0  0
  ##   I6T2  0  0  0  0  0 -1  0  0
  ##   I7T2  0  0  0  0  0  0 -1  0
  ##   I8T2  0  0  0  0  0  0  0 -1

# estimate model
# set intercept of second dimension (T2) to zero
beta.fixed <- cbind( 1, 2, 0 )
mod1c <- TAM::tam.mml( resp=dat, Q=Q, A=A1, beta.fixed=beta.fixed)
summary(mod1c)

#*************************************************
# Model 2: 2-dimensional 2PL model
#*************************************************

# set variance at T2 to 1
variance.fixed <- cbind(2,2,1)

# M2a: rasch.mml2 (in sirt)
mod2a <- sirt::rasch.mml2(dat, Q=Q, est.b=itemnr, est.a=itemnr, mu.fixed=mu.fixed,
             variance.fixed=variance.fixed, mmliter=100)
summary(mod2a)

#*************************************************
# Model 3: Concurrent calibration by assuming invariant item parameters
#*************************************************

library(mirt)   # use mirt for concurrent calibration
data(data.long)
dat <- data.long[,-1]
I <- ncol(dat)

# create user defined function for between item dimensionality 4PL model
name <- "4PLbw"
par <- c("low"=0,"upp"=1,"a"=1,"d"=0,"dimItem"=1)
est <- c(TRUE, TRUE,TRUE,TRUE,FALSE)
# item response function
irf <- function(par,Theta,ncat){
     low <- par[1]
     upp <- par[2]
     a <- par[3]
     d <- par[4]
     dimItem <- par[5]
     P1 <- low + ( upp - low ) * plogis( a*Theta[,dimItem] + d )
     cbind(1-P1, P1)
}

# create item response function
fourPLbetw <- mirt::createItem(name, par=par, est=est, P=irf)
head(dat)

# create mirt model (use variable names in mirt.model)
mirtsyn <- "
     T1=I1T1,I2T1,I3T1,I4T1,I5T1,I6T1
     T2=I3T2,I4T2,I5T2,I6T2,I7T2,I8T2
     COV=T1*T2,,T2*T2
     MEAN=T1
     CONSTRAIN=(I3T1,I3T2,d),(I4T1,I4T2,d),(I5T1,I5T2,d),(I6T1,I6T2,d),
                 (I3T1,I3T2,a),(I4T1,I4T2,a),(I5T1,I5T2,a),(I6T1,I6T2,a)
        "
# create mirt model
mirtmodel <- mirt::mirt.model( mirtsyn, itemnames=colnames(dat) )
# define parameters to be estimated
mod3.pars <- mirt::mirt(dat, mirtmodel$model, rep( "4PLbw",I),
                   customItems=list("4PLbw"=fourPLbetw), pars="values")
# select dimensions
ind <- intersect( grep("T2",mod3.pars$item), which( mod3.pars$name=="dimItem" ) )
mod3.pars[ind,"value"] <- 2
# set item parameters low and upp to non-estimated
ind <- which( mod3.pars$name \%in\% c("low","upp") )
mod3.pars[ind,"est"] <- FALSE

# estimate 2PL model
mod3 <- mirt::mirt(dat, mirtmodel$model, itemtype=rep( "4PLbw",I),
                customItems=list("4PLbw"=fourPLbetw), pars=mod3.pars, verbose=TRUE,
                technical=list(NCYCLES=50)  )
mirt.wrapper.coef(mod3)

#****** estimate model in lavaan
library(lavaan)

# specify syntax
lavmodel <- "
             #**** T1
             F1=~ a1*I1T1+a2*I2T1+a3*I3T1+a4*I4T1+a5*I5T1+a6*I6T1
             I1T1 | b1*t1 ; I2T1 | b2*t1 ; I3T1 | b3*t1 ; I4T1 | b4*t1
             I5T1 | b5*t1 ; I6T1 | b6*t1
             F1 ~~ 1*F1
             #**** T2
             F2=~ a3*I3T2+a4*I4T2+a5*I5T2+a6*I6T2+a7*I7T2+a8*I8T2
             I3T2 | b3*t1 ; I4T2 | b4*t1 ; I5T2 | b5*t1 ; I6T2 | b6*t1
             I7T2 | b7*t1 ; I8T2 | b8*t1
             F2 ~~ NA*F2
             F2 ~ 1
             #*** covariance
             F1 ~~ F2
                "
# estimate model using theta parameterization
mod3lav <- lavaan::cfa( data=dat, model=lavmodel,
            std.lv=TRUE, ordered=colnames(dat), parameterization="theta")
summary(mod3lav, standardized=TRUE, fit.measures=TRUE, rsquare=TRUE)

#*************************************************
# Model 4: Linking with items of different item slope groups
#*************************************************

data(data.long)
dat <- data.long
# dataset for T1
dat1 <- dat[, grep( "T1", colnames(dat) ) ]
colnames(dat1) <- gsub("T1","", colnames(dat1) )
# dataset for T2
dat2 <- dat[, grep( "T2", colnames(dat) ) ]
colnames(dat2) <- gsub("T2","", colnames(dat2) )

# 2PL model with slope groups T1
mod1 <- sirt::rasch.mml2( dat1, est.a=c( rep(1,2), rep(2,4) ) )
summary(mod1)

# 2PL model with slope groups T2
mod2 <- sirt::rasch.mml2( dat2, est.a=c( rep(1,4), rep(2,2) ) )
summary(mod2)

#------- Link 1: Haberman Linking
# collect item parameters
dfr1 <- data.frame( "study1", mod1$item$item, mod1$item$a, mod1$item$b )
dfr2 <- data.frame( "study2", mod2$item$item, mod2$item$a, mod2$item$b )
colnames(dfr2) <- colnames(dfr1) <- c("study", "item", "a", "b" )
itempars <- rbind( dfr1, dfr2 )
# Linking
link1 <- sirt::linking.haberman(itempars=itempars)

#------- Link 2: Invariance alignment method
# create objects for invariance.alignment
nu <- rbind( c(mod1$item$thresh,NA,NA), c(NA,NA,mod2$item$thresh) )
lambda <- rbind( c(mod1$item$a,NA,NA), c(NA,NA,mod2$item$a ) )
colnames(lambda) <- colnames(nu) <- paste0("I",1:8)
rownames(lambda) <- rownames(nu) <- c("T1", "T2")
# Linking
link2a <- sirt::invariance.alignment( lambda, nu )
summary(link2a)
}
}
%% \keyword{datasets}
