\name{data.big5}
\alias{data.big5}
\alias{data.big5.qgraph}


\docType{data}
\title{
Dataset Big 5 from \pkg{qgraph} Package
}

\description{
This is a Big 5 dataset from the \pkg{qgraph} package (Dolen, Oorts,
Stoel, Wicherts, 2009). It contains 500 subjects on 240 items.
}

\usage{
data(data.big5)
data(data.big5.qgraph)
}

\format{
\itemize{
\item The format of \code{data.big5} is: \cr 
\code{ num [1:500, 1:240] 1 0 0 0 0 1 1 2 0 1 ...} \cr
\code{ - attr(*, "dimnames")=List of 2} \cr
\code{  ..$ : NULL} \cr
\code{  ..$ : chr [1:240] "N1" "E2" "O3" "A4" ...} \cr

\item The format of \code{data.big5.qgraph} is: \cr  

\code{ num [1:500, 1:240] 2 3 4 4 5 2 2 1 4 2 ...} \cr
\code{ - attr(*, "dimnames")=List of 2} \cr
\code{  ..$ : NULL} \cr
\code{  ..$ : chr [1:240] "N1" "E2" "O3" "A4" ...} \cr
}
}


\details{
In these datasets, there exist 48 items for each dimension. The Big 5
dimensions are Neuroticism (\code{N}), Extraversion (\code{E}), 
Openness (\code{O}), Agreeableness (\code{A}) and 
Conscientiousness (\code{C}). Note that the \code{data.big5} differs from
\code{data.big5.qgraph} in a way that original items were recoded into
three categories 0,1 and 2.
}


\source{
See \code{big5} in \pkg{qgraph} package.
}


\references{
Dolan, C. V., Oort, F. J., Stoel, R. D., & Wicherts, J. M. (2009). 
Testing measurement invariance in the target rotates multigroup exploratory 
factor model. \emph{Structural Equation Modeling}, \bold{16}, 295-314.
}






\examples{
\dontrun{
# list of needed packages for the following examples
packages <- scan(what="character")
     sirt   TAM   eRm   CDM   mirt  ltm   mokken  psychotools  psychomix
     psych
   
# load packages. make an installation if necessary     
miceadds::library_install(packages)

#############################################################################
# EXAMPLE 1: Unidimensional models openness scale
#############################################################################

data(data.big5)	
# extract first 10 openness items
items <- which( substring( colnames(data.big5) , 1 , 1 ) == "O"  )[1:10]
dat <- data.big5[ , items ]
I <- ncol(dat)
summary(dat)
  ##   > colnames(dat)
  ##    [1] "O3"  "O8"  "O13" "O18" "O23" "O28" "O33" "O38" "O43" "O48"
# descriptive statistics
psych::describe(dat)
        
#****************
# Model 1: Partial credit model
#****************

#-- M1a: rm.facets (in sirt)
m1a <- sirt::rm.facets( dat )
summary(m1a)

#-- M1b: tam.mml (in TAM)
m1b <- TAM::tam.mml( resp=dat )
summary(m1b)

#-- M1c: gdm (in CDM)
theta.k <- seq(-6,6,len=21)
m1c <- CDM::gdm( dat , irtmodel="1PL" ,theta.k=theta.k , skillspace="normal")
summary(m1c)
# compare results with loglinear skillspace
m1c2 <- CDM::gdm( dat , irtmodel="1PL" ,theta.k=theta.k , skillspace="loglinear")
summary(m1c2)

#-- M1d: PCM (in eRm)
m1d <- eRm::PCM( dat )
summary(m1d)

#-- M1e: gpcm (in ltm)
m1e <- ltm::gpcm( dat , constraint = "1PL" , control=list(verbose=TRUE))
summary(m1e)

#-- M1f: mirt (in mirt)
m1f <- mirt::mirt( dat , model=1 , itemtype="1PL" , verbose=TRUE)
summary(m1f)
coef(m1f)

#-- M1g: PCModel.fit (in psychotools)
mod1g <- psychotools::PCModel.fit(dat)
summary(mod1g)
plot(mod1g)

#****************
# Model 2: Generalized partial credit model
#****************

#-- M2a: rm.facets (in sirt)
m2a <- sirt::rm.facets( dat , est.a.item=TRUE)
summary(m2a)
# Note that in rm.facets the mean of item discriminations is fixed to 1

#-- M2b: tam.mml.2pl (in TAM)
m2b <- TAM::tam.mml.2pl( resp=dat , irtmodel="GPCM")
summary(m2b)

#-- M2c: gdm (in CDM)
m2c <- CDM::gdm( dat , irtmodel="2PL" ,theta.k=seq(-6,6,len=21) ,
                   skillspace="normal" , standardized.latent=TRUE)
summary(m2c)

#-- M2d: gpcm (in ltm)
m2d <- ltm::gpcm( dat , control=list(verbose=TRUE))
summary(m2d)

#-- M2e: mirt (in mirt)
m2e <- mirt::mirt( dat , model=1 ,  itemtype="GPCM" , verbose=TRUE)
summary(m2e)
coef(m2e)

#****************
# Model 3: Nonparametric item response model
#****************

#-- M3a: ISOP and ADISOP model - isop.poly (in sirt)
m3a <- sirt::isop.poly( dat )
summary(m3a)
plot(m3a)

#-- M3b: Mokken scale analysis (in mokken) 
# Scalability coefficients
mokken::coefH(dat)
# Assumption of monotonicity
monotonicity.list <- mokken::check.monotonicity(dat)
summary(monotonicity.list)
plot(monotonicity.list)
# Assumption of non-intersecting ISRFs using method restscore
restscore.list <- mokken::check.restscore(dat)
summary(restscore.list)
plot(restscore.list) 

#****************
# Model 4: Graded response model
#****************

#-- M4a: mirt (in mirt)
m4a <- mirt::mirt( dat , model=1 ,  itemtype="graded" , verbose=TRUE)
print(m4a)
mirt.wrapper.coef(m4a)

#----  M4b: WLSMV estimation with cfa (in lavaan)
lavmodel <- "F =~ O3__O48
             F ~~ 1*F  
                "
# transform lavaan syntax with lavaanify.IRT
lavmodel <- TAM::lavaanify.IRT( lavmodel , items=colnames(dat) )$lavaan.syntax
mod4b <- lavaan::cfa( data= as.data.frame(dat) , model=lavmodel, std.lv = TRUE, 
                 ordered=colnames(dat) ,  parameterization="theta")
summary(mod4b , standardized=TRUE , fit.measures=TRUE , rsquare=TRUE)
coef(mod4b)

#****************
# Model 5: Normally distributed residuals
#****************

#----  M5a: cfa (in lavaan)
lavmodel <- "F =~ O3__O48
             F ~~ 1*F  
             F ~ 0*1
             O3__O48 ~ 1
                "                               
lavmodel <- TAM::lavaanify.IRT( lavmodel , items=colnames(dat) )$lavaan.syntax
mod5a <- lavaan::cfa( data= as.data.frame(dat) , model=lavmodel, std.lv = TRUE ,
                 estimator="MLR" )
summary(mod5a , standardized=TRUE , fit.measures=TRUE , rsquare=TRUE)

#----  M5b: mirt (in mirt)

# create user defined function
name <- 'normal'
par <- c("d" = 1 , "a1" = 0.8 , "vy" = 1)
est <- c(TRUE, TRUE,FALSE)
P.normal <- function(par,Theta,ncat){
     d <- par[1]
     a1 <- par[2]
     vy <- par[3]     
     psi <- vy - a1^2     
     # expected values given Theta
     mui <- a1*Theta[,1] + d
     TP <- nrow(Theta)
     probs <- matrix( NA , nrow=TP, ncol= ncat )
     eps <- .01
     for (cc in 1:ncat){
        probs[,cc] <- stats::dnorm( cc , mean = mui   , sd = sqrt( abs( psi + eps) ) )
                    }
     psum <- matrix( rep(rowSums( probs ),each=ncat)  , nrow=TP , ncol=ncat , byrow=TRUE)    
     probs <- probs / psum
     return(probs)                    
}

# create item response function
normal <- mirt::createItem(name, par=par, est=est, P=P.normal)
customItems <- list("normal"=normal) 
itemtype <- rep( "normal",I)
# define parameters to be estimated
mod5b.pars <- mirt::mirt(dat, 1, itemtype=itemtype , 
                   customItems=customItems , pars = "values")                
ind <- which( mod5b.pars$name == "vy")
vy <- apply( dat , 2 , var , na.rm=TRUE ) 
mod5b.pars[ ind , "value" ] <- vy
ind <- which( mod5b.pars$name == "a1")
mod5b.pars[ ind , "value" ] <- .5* sqrt(vy)
ind <- which( mod5b.pars$name == "d")
mod5b.pars[ ind , "value" ] <- colMeans( dat , na.rm=TRUE )
                                                       
# estimate model
mod5b <- mirt::mirt(dat, 1, itemtype=itemtype , customItems=customItems ,  
                 pars = mod5b.pars , verbose=TRUE    )
sirt::mirt.wrapper.coef(mod5b)$coef

# some item plots
    par(ask=TRUE)
plot(mod5b, type = 'trace', layout = c(1,1))
    par(ask=FALSE)
# Alternatively:
sirt::mirt.wrapper.itemplot(mod5b) 
}
}


\keyword{datasets}
