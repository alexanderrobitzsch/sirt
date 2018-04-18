%% File Name: data.mcdonald.Rd
%% File Version: 0.15

\name{data.mcdonald}

\alias{data.mcdonald.act15}
\alias{data.mcdonald.LSAT6}
\alias{data.mcdonald.rape}

\docType{data}

\title{
Some Datasets from McDonald's \emph{Test Theory} Book
}

\description{
Some datasets from McDonald (1999), especially related to using
NOHARM for item response modeling. See Examples below.
}


\usage{
data(data.mcdonald.act15)
data(data.mcdonald.LSAT6)
data(data.mcdonald.rape)
}



\format{
\itemize{	
\item The format of the ACT15 data \code{data.mcdonald.act15} is: 

\code{ num [1:15, 1:15] 0.49 0.44 0.38 0.3 0.29 0.13 0.23 0.16 0.16 0.23 ...} \cr
\code{ - attr(*, "dimnames")=List of 2} \cr
\code{  ..$ : chr [1:15] "A01" "A02" "A03" "A04" ...} \cr
\code{  ..$ : chr [1:15] "A01" "A02" "A03" "A04" ...}

The dataset (which is the product-moment covariance matrix)
is obtained from Ch. 12 in McDonald (1999). \cr

\item The format of the LSAT6 data \code{data.mcdonald.LSAT6} is: 

\code{'data.frame':   1004 obs. of  5 variables:} \cr
\code{ $ L1: int  0 0 0 0 0 0 0 0 0 0 ...} \cr
\code{ $ L2: int  0 0 0 0 0 0 0 0 0 0 ...} \cr
\code{ $ L3: int  0 0 0 0 0 0 0 0 0 0 ...} \cr
\code{ $ L4: int  0 0 0 0 0 0 0 0 0 1 ...} \cr
\code{ $ L5: int  0 0 0 1 1 1 1 1 1 0 ...}

The dataset is obtained from Ch. 6 in McDonald (1999). \cr

\item The format of the rape myth scale data  \code{data.mcdonald.rape} is

\code{List of 2} \cr
\code{ $ lambda: num [1:2, 1:19] 1.13 0.88 0.85 0.77 0.79 0.55 1.12 1.01 0.99 0.79 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 2} \cr
\code{  .. ..$ : chr [1:2] "male" "female"} \cr
\code{  .. ..$ : chr [1:19] "I1" "I2" "I3" "I4" ...} \cr
\code{ $ nu    : num [1:2, 1:19] 2.88 1.87 3.12 2.32 2.13 1.43 3.79 2.6 3.01 2.11 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 2} \cr
\code{  .. ..$ : chr [1:2] "male" "female"} \cr
\code{  .. ..$ : chr [1:19] "I1" "I2" "I3" "I4" ...} \cr

The dataset is obtained from Ch. 15 in McDonald (1999).

}
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Tables in McDonald (1999)
}


\references{
McDonald, R. P. (1999). \emph{Test theory: A unified treatment}. 
Psychology Press.
}


\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: LSAT6 data	| Chapter 12 McDonald (1999)
#############################################################################	
data(data.mcdonald.act15)

#************
# Model 1: 2-parameter normal ogive model

#++ NOHARM estimation
I <- ncol(dat)
# covariance structure
P.pattern <- matrix( 0 , ncol=1 , nrow=1 )
P.init <- 1+0*P.pattern
# fix all entries in the loading matrix to 1
F.pattern <- matrix( 1 , nrow=I , ncol=1 )
F.init <- F.pattern        
# estimate model
mod1a <- sirt::R2noharm( dat = dat , model.type="CFA" , F.pattern = F.pattern , 
             F.init = F.init , P.pattern = P.pattern , P.init = P.init ,
             writename = "LSAT6__1dim_2pno" , noharm.path = noharm.path , dec ="," )
summary(mod1a , logfile="LSAT6__1dim_2pno__SUMMARY")

#++ pairwise marginal maximum likelihood estimation using the probit link
mod1b <- sirt::rasch.pml3( dat , est.a=1:I , est.sigma=FALSE)

#************
# Model 2: 1-parameter normal ogive model

#++ NOHARM estimation
# covariance structure
P.pattern <- matrix( 0 , ncol=1 , nrow=1 )
P.init <- 1+0*P.pattern
# fix all entries in the loading matrix to 1
F.pattern <- matrix( 2 , nrow=I , ncol=1 )
F.init <- 1+0*F.pattern       
# estimate model
mod2a <- sirt::R2noharm( dat = dat , model.type="CFA" , F.pattern = F.pattern ,
                F.init = F.init , P.pattern = P.pattern , P.init = P.init , 
                writename = "LSAT6__1dim_1pno" , noharm.path = noharm.path , dec ="," )
summary(mod2a , logfile="LSAT6__1dim_1pno__SUMMARY")

# PMML estimation
mod2b <- sirt::rasch.pml3( dat , est.a=rep(1,I) , est.sigma=FALSE )
summary(mod2b)

#************
# Model 3: 3-parameter normal ogive model with fixed guessing parameters

#++ NOHARM estimation
# covariance structure
P.pattern <- matrix( 0 , ncol=1 , nrow=1 )
P.init <- 1+0*P.pattern
# fix all entries in the loading matrix to 1
F.pattern <- matrix( 1 , nrow=I , ncol=1 )
F.init <- 1+0*F.pattern      
# estimate model
mod <- sirt::R2noharm( dat = dat , model.type="CFA" ,  guesses=rep(.2,I) , 
            F.pattern = F.pattern , F.init = F.init , P.pattern = P.pattern ,
            P.init = P.init , writename = "LSAT6__1dim_3pno", 
            noharm.path = noharm.path , dec ="," )
summary(mod , logfile="LSAT6__1dim_3pno__SUMMARY")

#++ logistic link function employed in smirt function 
mod1d <- sirt::smirt(dat, Qmatrix=F.pattern, est.a= matrix(1:I,I,1), c.init=rep(.2,I))
summary(mod1d)

#############################################################################
# EXAMPLE 2: ACT15 data	| Chapter 6 McDonald (1999)
#############################################################################	
data(data.mcdonald.act15)
pm <- data.mcdonald.act15

#************
# Model 1: 2-dimensional exploratory factor analysis
mod1 <- sirt::R2noharm( pm=pm , n=1000, model.type="EFA" , dimensions=2 , 
             writename = "ACT15__efa_2dim" , noharm.path = noharm.path , dec ="," )
summary(mod1)

#************
# Model 2: 2-dimensional independent clusters basis solution
P.pattern <- matrix(1,2,2)
diag(P.pattern) <- 0
P.init <- 1+0*P.pattern
F.pattern <- matrix(0,15,2)
F.pattern[ c(1:5,11:15),1] <- 1
F.pattern[ c(6:10,11:15),2] <- 1
F.init <- F.pattern

# estimate model
mod2 <- sirt::R2noharm( pm=pm , n=1000 ,  model.type="CFA" , F.pattern = F.pattern , 
            F.init = F.init , P.pattern = P.pattern ,P.init = P.init , 
            writename = "ACT15_indep_clusters" , noharm.path = noharm.path , dec ="," )
summary(mod2)

#************
# Model 3: Hierarchical model

P.pattern <- matrix(0,3,3)
P.init <- P.pattern
diag(P.init) <- 1
F.pattern <- matrix(0,15,3)
F.pattern[,1] <- 1	# all items load on g factor
F.pattern[ c(1:5,11:15),2] <- 1   # Items 1-5 and 11-15 load on first nested factor
F.pattern[ c(6:10,11:15),3] <- 1  # Items 6-10 and 11-15 load on second nested factor
F.init <- F.pattern

# estimate model
mod3 <- sirt::R2noharm( pm=pm , n=1000 ,  model.type="CFA" , F.pattern = F.pattern , 
           F.init = F.init , P.pattern = P.pattern , P.init = P.init , 
           writename = "ACT15_hierarch_model" , noharm.path = noharm.path , dec ="," )
summary(mod3)

#############################################################################
# EXAMPLE 3: Rape myth scale | Chapter 15 McDonald (1999)
#############################################################################   
data(data.mcdonald.rape)
lambda <- data.mcdonald.rape$lambda
nu <- data.mcdonald.rape$nu

# obtain multiplier for factor loadings (Formula 15.5)
k <- sum( lambda[1,] * lambda[2,] ) / sum( lambda[2,]^2 )
  ##   [1] 1.263243

# additive parameter (Formula 15.7)
c <- sum( lambda[2,]*(nu[1,]-nu[2,]) ) / sum( lambda[2,]^2 )
  ##   [1] 1.247697

# SD in the female group
1/k
  ##   [1] 0.7916132

# M in the female group
- c/k
  ##   [1] -0.9876932

# Burt's coefficient of factorial congruence (Formula 15.10a)
sum( lambda[1,] * lambda[2,] ) / sqrt( sum( lambda[1,]^2 ) * sum( lambda[2,]^2 ) )
  ##   [1] 0.9727831

# congruence for mean parameters
sum(  (nu[1,]-nu[2,]) * lambda[2,] ) / sqrt( sum( (nu[1,]-nu[2,])^2 ) * sum( lambda[2,]^2 ) )
  ##   [1] 0.968176
}
}
\keyword{datasets}
