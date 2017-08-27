%% File Name: nedelsky.sim.Rd
%% File Version: 0.16
%% File Last Change: 2017-01-18 18:19:11

\name{nedelsky-methods}

\alias{nedelsky.sim}
\alias{nedelsky.latresp}
\alias{nedelsky.irf}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions for the Nedelsky Model
}

\description{
Functions for simulating and estimating the Nedelsky model
(Bechger et al., 2003, 2005). \code{nedelsky.sim} can be used for
simulating the model, \code{nedelsky.irf} computes the item response
function and can be used for example when estimating the
Nedelsky model in the \pkg{mirt} package.
}

\usage{
# simulating the Nedelsky model	
nedelsky.sim(theta, b, a = NULL, tau = NULL)

# creating latent responses of the Nedelsky model
nedelsky.latresp(K)

# computing the item response function of the Nedelsky model
nedelsky.irf(Theta, K, b, a, tau, combis, thdim = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
Unidimensional ability (theta)
}
  \item{b}{
Matrix of category difficulties
}
  \item{a}{
Vector of item discriminations
}
  \item{tau}{
Category attractivity parameters \eqn{\tau} (see Bechger et al., 2005) 
}
\item{K}{(Maximum) Number of distractors of the used multiple choice items}
  \item{Theta}{
Theta vector. Note that the Nedelsky model can be only specified
  as models with between item dimensionality (defined in \code{thdim}).
}
  \item{combis}{
Latent response classes as produced by \code{nedelsky.latresp}.
}
  \item{thdim}{
Theta dimension at which the item loads
}
}

\details{
Assume that for item \eqn{i} there exists \eqn{K+1} categories \eqn{0,1,...,K}.
The category 0 denotes the correct alternative. The Nedelsky model assumes
that a respondent eliminates all distractors which are thought to be
incorrect and guesses the solution from the remaining alternatives.
This means, that for item \eqn{i}, \eqn{K} latent variables \eqn{S_{ik}}
are defined which indicate whether alternative \eqn{k} has been correctly
identified as a distractor. By definition, the correct alternative is never
been judged as wrong by the respondent.

Formally, the Nedelsky model assumes a 2PL model for eliminating each of the 
distractors
\deqn{P(S_{ik}=1 | \theta ) = invlogit[ a_i ( \theta - b_{ik} ) ] }
where \eqn{\theta} is the person ability and \eqn{b_{ik}} are
distractor difficulties.

The guessing process of the Nedelsky model is defined as
\deqn{P(X_i= j | \theta , S_{i1} , ... , S_{iK} ) =
\frac{ ( 1- S_{ij} ) \tau_{ij} }{ \sum_{k=0}^K [ ( 1- S_{ik} ) \tau_{ik} ] }}
where \eqn{\tau_{ij}} are attractivity parameters of alternative \eqn{j}.
By definition \eqn{\tau_{i0}} is set to 1. By default, all attractivity parameters
are set to 1.
}



%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}


\references{
Bechger, T. M., Maris, G., Verstralen, H. H. F. M., & Verhelst, N. D. (2003). 
\emph{The Nedelsky model for multiple-choice items}. 
CITO Research Report, 2003-5.

Bechger, T. M., Maris, G., Verstralen, H. H. F. M., & Verhelst, N. D. (2005). 
The Nedelsky model for multiple-choice items. 
In L. van der Ark, M. Croon, & Sijtsma, K. (Eds.).
\emph{New developments in categorical data analysis for the social and behavioral 
sciences}, pp. 187-206. Mahwah, Lawrence Erlbaum. 
}


\author{
Alexander Robitzsch
}

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
# EXAMPLE 1: Simulated data according to the Nedelsky model
#############################################################################

#*** simulate data
set.seed(123)
I <- 20          # number of items
b <- matrix(NA,I,ncol=3)
b[,1] <- -0.5 + stats::runif( I , -.75 , .75 )
b[,2] <- -1.5 + stats::runif( I , -.75 , .75 )
b[,3] <- -2.5 + stats::runif( I , -.75 , .75 )
K <- 3           # number of distractors
N <- 2000        # number of persons
# apply simulation function
dat <- sirt::nedelsky.sim( theta = stats::rnorm(N,sd=1.2) , b=b )

#*** latent response patterns
K <- 3
combis <- sirt::nedelsky.latresp(K=3)

#*** defining the Nedelsky item response function for estimation in mirt
par <- c( 3  , rep(-1,K) , 1 , rep(1,K+1) ,1)
names(par) <- c("K" , paste0("b",1:K) , "a" , paste0("tau" , 0:K) ,"thdim")
est <- c( FALSE , rep(TRUE,K) , rep(FALSE , K+1 + 2 ) )
names(est) <- names(par)
nedelsky.icc <- function( par , Theta , ncat ){
     K <- par[1]
     b <- par[ 1:K + 1]
     a <- par[ K+2]
     tau <- par[1:(K+1) + (K+2) ]
     thdim <- par[ K+2+K+1 +1 ]
     probs <- sirt::nedelsky.irf( Theta , K=K , b=b , a=a , tau=tau , combis , 
                    thdim=thdim  )$probs
     return(probs)
}
name <- "nedelsky"
# create item response function
nedelsky.itemfct <- mirt::createItem(name, par=par, est=est, P=nedelsky.icc)

#*** define model in mirt
mirtmodel <- mirt::mirt.model("
           F1 = 1-20
           COV = F1*F1
           # define some prior distributions
           PRIOR = (1-20,b1 ,norm,-1,2),(1-20,b2 ,norm,-1,2),
                   (1-20,b3 ,norm,-1,2)
        " )

itemtype <- rep("nedelsky" , I )
customItems <- list("nedelsky"= nedelsky.itemfct)
# define parameters to be estimated
mod1.pars <- mirt::mirt(dat, mirtmodel , itemtype=itemtype , 
                   customItems=customItems, pars = "values")                
# estimate model
mod1 <- mirt::mirt(dat,mirtmodel , itemtype=itemtype , customItems=customItems, 
               pars = mod1.pars , verbose=TRUE  )
# model summaries
print(mod1)                
summary(mod1)
mirt.wrapper.coef( mod1 )$coef
mirt.wrapper.itemplot(mod1 ,ask=TRUE)

#******************************************************
# fit Nedelsky model with xxirt function in sirt

# define item class for xxirt
item_nedelsky <- sirt::xxirt_createDiscItem( name = "nedelsky" , par = par , 
                est = est , P = nedelsky.icc ,
                prior = c( b1="dnorm" , b2 = "dnorm" , b3 = "dnorm" ) , 
                prior_par1 = c( b1 = -1 , b2=-1 , b3=-1) , 
                prior_par2 = c(b1=2, b2=2 , b3=2) )
customItems <- list( item_nedelsky )

#---- definition theta distribution
#** theta grid
Theta <- matrix( seq(-6,6,length=21) , ncol=1 )
#** theta distribution
P_Theta1 <- function( par , Theta , G){
    mu <- par[1]
    sigma <- max( par[2] , .01 )
    TP <- nrow(Theta)
    pi_Theta <- matrix( 0 , nrow=TP , ncol=G)
    pi1 <- dnorm( Theta[,1] , mean = mu , sd = sigma )
    pi1 <- pi1 / sum(pi1)
    pi_Theta[,1] <- pi1
    return(pi_Theta)
}
#** create distribution class
par_Theta <- c( "mu"=0, "sigma" = 1 )  
customTheta <- sirt::xxirt_createThetaDistribution( par=par_Theta, est=c(FALSE,TRUE), 
                   P=P_Theta1 )

#-- create parameter table
itemtype <- rep( "nedelsky" , I )
partable <- sirt::xxirt_createParTable( dat, itemtype = itemtype, customItems = customItems)

# estimate model
mod2 <- sirt::xxirt( dat=dat, Theta=Theta, partable=partable, customItems=customItems,
                    customTheta = customTheta)
summary(mod2)
# compare sirt::xxirt and mirt::mirt
logLik(mod2)
mod1@Fit$logLik

#############################################################################
# EXAMPLE 2: Multiple choice dataset data.si06
#############################################################################

data(data.si06)
dat <- data.si06

#*** create latent responses
combis <- sirt::nedelsky.latresp(K)
I <- ncol(dat)
#*** define item response function
K <- 3
par <- c( 3  , rep(-1,K) , 1 , rep(1,K+1) ,1)
names(par) <- c("K" , paste0("b",1:K) , "a" , paste0("tau" , 0:K) ,"thdim")
est <- c( FALSE , rep(TRUE,K) , rep(FALSE , K+1 + 2 ) )
names(est) <- names(par)
nedelsky.icc <- function( par , Theta , ncat ){
     K <- par[1]
     b <- par[ 1:K + 1]
     a <- par[ K+2]
     tau <- par[1:(K+1) + (K+2) ]
     thdim <- par[ K+2+K+1 +1 ]
     probs <- sirt::nedelsky.irf( Theta , K=K , b=b , a=a , tau=tau , combis , 
                    thdim=thdim  )$probs
     return(probs)
}
name <- "nedelsky"
# create item response function
nedelsky.itemfct <- mirt::createItem(name, par=par, est=est, P=nedelsky.icc)

#*** define model in mirt
mirtmodel <- mirt::mirt.model("
           F1 = 1-14
           COV = F1*F1
           PRIOR = (1-14,b1 ,norm,-1,2),(1-14,b2 ,norm,-1,2),
                   (1-14,b3 ,norm,-1,2)
        " )

itemtype <- rep("nedelsky" , I )
customItems <- list("nedelsky"= nedelsky.itemfct)
# define parameters to be estimated
mod1.pars <- mirt::mirt(dat, mirtmodel , itemtype=itemtype , 
                   customItems=customItems, pars = "values")                

#*** estimate model
mod1 <- mirt::mirt(dat,mirtmodel , itemtype=itemtype , customItems=customItems, 
               pars = mod1.pars , verbose=TRUE )
#*** summaries                
print(mod1)                
summary(mod1)
mirt.wrapper.coef( mod1 )$coef
mirt.wrapper.itemplot(mod1 ,ask=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{Nedelsky model}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
