%% File Name: data.read.Rd
%% File Version: 1.997

\name{data.read}
\alias{data.read}
\docType{data}
\title{
Dataset Reading
}
\description{
This dataset contains \eqn{N=328} students and \eqn{I=12} items measuring reading
competence. All 12 items are arranged into 3 testlets (items with common
text stimulus) labeled as
A, B and C. The allocation of items to testlets is indicated by their
variable names.
}
\usage{data(data.read)}
\format{
  A data frame with 328 persons on the following 12 variables.
  Rows correspond to persons and columns to items. The following items are
  included in \code{data.read}:

  Testlet A: \code{A1}, \code{A2}, \code{A3}, \code{A4}

  Testlet B: \code{B1}, \code{B2}, \code{B3}, \code{B4}

  Testlet C: \code{C1}, \code{C2}, \code{C3}, \code{C4}

}



\examples{
\dontrun{
data(data.read)
dat <- data.read
I <- ncol(dat)

# list of needed packages for the following examples
packages <- scan(what="character")
     eRm  ltm  TAM mRm  CDM  mirt psychotools  IsingFit  igraph  qgraph  pcalg
     poLCA  randomLCA psychomix MplusAutomation lavaan

# load packages. make an installation if necessary
miceadds::library_install(packages)

#*****************************************************
# Model 1: Rasch model
#*****************************************************

#----  M1a: rasch.mml2 (in sirt)
mod1a <- sirt::rasch.mml2(dat)
summary(mod1a)

#----  M1b: smirt (in sirt)
Qmatrix <- matrix(1,nrow=I, ncol=1)
mod1b <- sirt::smirt(dat,Qmatrix=Qmatrix)
summary(mod1b)

#----  M1c: gdm (in CDM)
theta.k <- seq(-6,6,len=21)
mod1c <- CDM::gdm(dat,theta.k=theta.k,irtmodel="1PL", skillspace="normal")
summary(mod1c)

#----  M1d: tam.mml (in TAM)
mod1d <- TAM::tam.mml( resp=dat )
summary(mod1d)

#----  M1e: RM (in eRm)
mod1e <- eRm::RM( dat )
  # eRm uses Conditional Maximum Likelihood (CML) as the estimation method.
summary(mod1e)
eRm::plotPImap(mod1e)

#----  M1f: mrm (in mRm)
mod1f <- mRm::mrm( dat, cl=1)   # CML estimation
mod1f$beta  # item parameters

#----  M1g: mirt (in mirt)
mod1g <- mirt::mirt( dat, model=1, itemtype="Rasch", verbose=TRUE )
print(mod1g)
summary(mod1g)
coef(mod1g)
    # arrange coefficients in nicer layout
sirt::mirt.wrapper.coef(mod1g)$coef

#----  M1h: rasch (in ltm)
mod1h <- ltm::rasch( dat, control=list(verbose=TRUE ) )
summary(mod1h)
coef(mod1h)

#----  M1i: RaschModel.fit (in psychotools)
mod1i <- psychotools::RaschModel.fit(dat)  # CML estimation
summary(mod1i)
plot(mod1i)

#----  M1j: noharm.sirt (in sirt)
Fpatt <- matrix( 0, I, 1 )
Fval <- 1 + 0*Fpatt
Ppatt <- Pval <- matrix(1,1,1)
mod1j <- sirt::noharm.sirt( dat=dat, Ppatt=Ppatt, Fpatt=Fpatt, Fval=Fval, Pval=Pval)
summary(mod1j)
  #   Normal-ogive model, multiply item discriminations with constant D=1.7.
  #   The same holds for other examples with noharm.sirt and R2noharm.
plot(mod1j)

#----  M1k: rasch.pml3 (in sirt)
mod1k <- sirt::rasch.pml3( dat=dat)
  #         pairwise marginal maximum likelihood estimation
summary(mod1k)

#----  M1l: running Mplus (using MplusAutomation package)
mplus_path <- "c:/Mplus7/Mplus.exe"    # locate Mplus executable
#****************
  # specify Mplus object
mplusmod <- MplusAutomation::mplusObject(
    TITLE="1PL in Mplus ;",
    VARIABLE=paste0( "CATEGORICAL ARE ", paste0(colnames(dat),collapse=" ") ),
    MODEL="
       ! fix all item loadings to 1
       F1 BY A1@1 A2@1 A3@1 A4@1 ;
       F1 BY B1@1 B2@1 B3@1 B4@1 ;
       F1 BY C1@1 C2@1 C3@1 C4@1 ;
       ! estimate variance
       F1 ;
            ",
    ANALYSIS="ESTIMATOR=MLR;",
    OUTPUT="stand;",
    usevariables=colnames(dat),  rdata=dat )
#****************

  # write Mplus syntax
filename <- "mod1u"   # specify file name
  # create Mplus syntaxes
res2 <- MplusAutomation::mplusModeler(object=mplusmod, dataout=paste0(filename,".dat"),
               modelout=paste0(filename,".inp"), run=0 )
  # run Mplus model
MplusAutomation::runModels( filefilter=paste0(filename,".inp"), Mplus_command=mplus_path)
  # alternatively, the system() command can also be used
  # get results
mod1l <- MplusAutomation::readModels(target=getwd(), filefilter=filename )
mod1l$summaries    # summaries
mod1l$parameters$unstandardized   # parameter estimates

#*****************************************************
# Model 2: 2PL model
#*****************************************************

#----  M2a: rasch.mml2 (in sirt)
mod2a <- sirt::rasch.mml2(dat, est.a=1:I)
summary(mod2a)

#----  M2b: smirt (in sirt)
mod2b <- sirt::smirt(dat,Qmatrix=Qmatrix,est.a="2PL")
summary(mod2b)

#----  M2c: gdm (in CDM)
mod2c <- CDM::gdm(dat,theta.k=theta.k,irtmodel="2PL", skillspace="normal")
summary(mod2c)

#----  M2d: tam.mml (in TAM)
mod2d <- TAM::tam.mml.2pl( resp=dat )
summary(mod2d)

#----  M2e: mirt (in mirt)
mod2e <- mirt::mirt( dat, model=1, itemtype="2PL" )
print(mod2e)
summary(mod2e)
sirt::mirt.wrapper.coef(mod1g)$coef

#----  M2f: ltm (in ltm)
mod2f <- ltm::ltm( dat ~ z1, control=list(verbose=TRUE ) )
summary(mod2f)
coef(mod2f)
plot(mod2f)

#----  M2g: R2noharm (in NOHARM, running from within R using sirt package)
  # define noharm.path where 'NoharmCL.exe' is located
noharm.path <- "c:/NOHARM"
  # covariance matrix
P.pattern <- matrix( 1, ncol=1, nrow=1 )
P.init <- P.pattern
P.init[1,1] <- 1
  # loading matrix
F.pattern <- matrix(1,I,1)
F.init <- F.pattern
  # estimate model
mod2g <- sirt::R2noharm( dat=dat, model.type="CFA", F.pattern=F.pattern,
             F.init=F.init, P.pattern=P.pattern, P.init=P.init,
             writename="ex2g", noharm.path=noharm.path, dec="," )
summary(mod2g)

#----  M2h: noharm.sirt (in sirt)
mod2h <- sirt::noharm.sirt( dat=dat, Ppatt=P.pattern,Fpatt=F.pattern,
              Fval=F.init, Pval=P.init )
summary(mod2h)
plot(mod2h)

#----  M2i: rasch.pml2 (in sirt)
mod2i <- sirt::rasch.pml2(dat, est.a=1:I)
summary(mod2i)

#----  M2j: WLSMV estimation with cfa (in lavaan)
lavmodel <- "F=~ A1+A2+A3+A4+B1+B2+B3+B4+
                        C1+C2+C3+C4"
mod2j <- lavaan::cfa( data=dat, model=lavmodel, std.lv=TRUE, ordered=colnames(dat))
summary(mod2j, standardized=TRUE, fit.measures=TRUE, rsquare=TRUE)

#*****************************************************
# Model 3: 3PL model (note that results can be quite unstable!)
#*****************************************************

#----  M3a: rasch.mml2 (in sirt)
mod3a <- sirt::rasch.mml2(dat, est.a=1:I, est.c=1:I)
summary(mod3a)

#----  M3b: smirt (in sirt)
mod3b <- sirt::smirt(dat,Qmatrix=Qmatrix,est.a="2PL", est.c=1:I)
summary(mod3b)

#----  M3c: mirt (in mirt)
mod3c <- mirt::mirt( dat, model=1, itemtype="3PL", verbose=TRUE)
summary(mod3c)
coef(mod3c)
  # stabilize parameter estimating using informative priors for guessing parameters
mirtmodel <- mirt::mirt.model("
            F=1-12
            PRIOR=(1-12, g, norm, -1.38, 0.25)
            ")
  # a prior N(-1.38,.25) is specified for transformed guessing parameters: qlogis(g)
  # simulate values from this prior for illustration
N <- 100000
logit.g <- stats::rnorm(N, mean=-1.38, sd=sqrt(.5) )
graphics::plot( stats::density(logit.g) )  # transformed qlogis(g)
graphics::plot( stats::density( stats::plogis(logit.g)) )  # g parameters
  # estimate 3PL with priors
mod3c1 <- mirt::mirt(dat, mirtmodel, itemtype="3PL",verbose=TRUE)
coef(mod3c1)
  # In addition, set upper bounds for g parameters of .35
mirt.pars <- mirt::mirt( dat, mirtmodel, itemtype="3PL",  pars="values")
ind <- which( mirt.pars$name=="g" )
mirt.pars[ ind, "value" ] <- stats::plogis(-1.38)
mirt.pars[ ind, "ubound" ] <- .35
  # prior distribution for slopes
ind <- which( mirt.pars$name=="a1" )
mirt.pars[ ind, "prior_1" ] <- 1.3
mirt.pars[ ind, "prior_2" ] <- 2
mod3c2 <- mirt::mirt(dat, mirtmodel, itemtype="3PL",
                pars=mirt.pars,verbose=TRUE, technical=list(NCYCLES=100) )
coef(mod3c2)
sirt::mirt.wrapper.coef(mod3c2)

#----  M3d: ltm (in ltm)
mod3d <- ltm::tpm( dat, control=list(verbose=TRUE), max.guessing=.3)
summary(mod3d)
coef(mod3d) #=> numerical instabilities

#*****************************************************
# Model 4: 3-dimensional Rasch model
#*****************************************************

# define Q-matrix
Q <- matrix( 0, nrow=12, ncol=3 )
Q[ cbind(1:12, rep(1:3,each=4) ) ] <- 1
rownames(Q) <- colnames(dat)
colnames(Q) <- c("A","B","C")

# define nodes
theta.k <- seq(-6,6,len=13)

#----  M4a: smirt (in sirt)
mod4a <- sirt::smirt(dat,Qmatrix=Q,irtmodel="comp", theta.k=theta.k, maxiter=30)
summary(mod4a)

#----  M4b: rasch.mml2 (in sirt)
mod4b <- sirt::rasch.mml2(dat,Q=Q,theta.k=theta.k, mmliter=30)
summary(mod4b)

#----  M4c: gdm (in CDM)
mod4c <- CDM::gdm( dat, irtmodel="1PL", theta.k=theta.k, skillspace="normal",
            Qmatrix=Q, maxiter=30, centered.latent=TRUE )
summary(mod4c)

#----  M4d: tam.mml (in TAM)
mod4d <- TAM::tam.mml( resp=dat, Q=Q, control=list(nodes=theta.k, maxiter=30) )
summary(mod4d)

#----  M4e: R2noharm (in NOHARM, running from within R using sirt package)
noharm.path <- "c:/NOHARM"
  # covariance matrix
P.pattern <- matrix( 1, ncol=3, nrow=3 )
P.init <- 0.8+0*P.pattern
diag(P.init) <- 1
  # loading matrix
F.pattern <- 0*Q
F.init <- Q
  # estimate model
mod4e <- sirt::R2noharm( dat=dat, model.type="CFA", F.pattern=F.pattern,
    F.init=F.init, P.pattern=P.pattern, P.init=P.init,
    writename="ex4e", noharm.path=noharm.path, dec="," )
summary(mod4e)

#----  M4f: mirt (in mirt)
cmodel <- mirt::mirt.model("
     F1=1-4
     F2=5-8
     F3=9-12
     # equal item slopes correspond to the Rasch model
     CONSTRAIN=(1-4, a1), (5-8, a2), (9-12,a3)
     COV=F1*F2, F1*F3, F2*F3
     " )
mod4f <- mirt::mirt(dat, cmodel, verbose=TRUE)
summary(mod4f)

#*****************************************************
# Model 5: 3-dimensional 2PL model
#*****************************************************

#----  M5a: smirt (in sirt)
mod5a <- sirt::smirt(dat,Qmatrix=Q,irtmodel="comp", est.a="2PL", theta.k=theta.k,
                 maxiter=30)
summary(mod5a)

#----  M5b: rasch.mml2 (in sirt)
mod5b <- sirt::rasch.mml2(dat,Q=Q,theta.k=theta.k,est.a=1:12, mmliter=30)
summary(mod5b)

#----  M5c: gdm (in CDM)
mod5c <- CDM::gdm( dat, irtmodel="2PL", theta.k=theta.k, skillspace="loglinear",
            Qmatrix=Q, maxiter=30, centered.latent=TRUE,
            standardized.latent=TRUE)
summary(mod5c)

#----  M5d: tam.mml (in TAM)
mod5d <- TAM::tam.mml.2pl( resp=dat, Q=Q, control=list(nodes=theta.k, maxiter=30) )
summary(mod5d)

#----  M5e: R2noharm (in NOHARM, running from within R using sirt package)
noharm.path <- "c:/NOHARM"
  # covariance matrix
P.pattern <- matrix( 1, ncol=3, nrow=3 )
diag(P.pattern) <- 0
P.init <- 0.8+0*P.pattern
diag(P.init) <- 1
  # loading matrix
F.pattern <- Q
F.init <- Q
  # estimate model
mod5e <- sirt::R2noharm( dat=dat, model.type="CFA", F.pattern=F.pattern,
    F.init=F.init, P.pattern=P.pattern, P.init=P.init,
    writename="ex5e", noharm.path=noharm.path, dec="," )
summary(mod5e)

#----  M5f: mirt (in mirt)
cmodel <- mirt::mirt.model("
   F1=1-4
   F2=5-8
   F3=9-12
   COV=F1*F2, F1*F3, F2*F3
   "  )
mod5f <- mirt::mirt(dat, cmodel, verbose=TRUE)
summary(mod5f)

#*****************************************************
# Model 6: Network models (Graphical models)
#*****************************************************

#----  M6a: Ising model using the IsingFit package (undirected graph)
#        - fit Ising model using the "OR rule" (AND=FALSE)
mod6a <- IsingFit::IsingFit(x=dat, family="binomial", AND=FALSE)
summary(mod6a)
##           Network Density:                 0.29
##    Gamma:                  0.25
##    Rule used:              Or-rule
# plot results
qgraph::qgraph(mod6a$weiadj,fade=FALSE)

#**-- graph estimation using pcalg package

# some packages from Bioconductor must be downloaded at first (if not yet done)
if (FALSE){  # set 'if (TRUE)' if packages should be downloaded
     source("http://bioconductor.org/biocLite.R")
     biocLite("RBGL")
     biocLite("Rgraphviz")
}

#----  M6b: graph estimation based on Pearson correlations
V <- colnames(dat)
n <- nrow(dat)
mod6b <- pcalg::pc(suffStat=list(C=stats::cor(dat), n=n ),
             indepTest=gaussCItest, ## indep.test: partial correlations
             alpha=0.05, labels=V, verbose=TRUE)
plot(mod6b)
# plot in qgraph package
qgraph::qgraph(mod6b, label.color=rep( c( "red", "blue","darkgreen" ), each=4 ),
         edge.color="black")
summary(mod6b)

#----  M6c: graph estimation based on tetrachoric correlations
mod6c <- pcalg::pc(suffStat=list(C=sirt::tetrachoric2(dat)$rho, n=n ),
             indepTest=gaussCItest, alpha=0.05, labels=V, verbose=TRUE)
plot(mod6c)
summary(mod6c)

#----  M6d: Statistical implicative analysis (in sirt)
mod6d <- sirt::sia.sirt(dat, significance=.85 )
  # plot results with igraph and qgraph package
plot( mod6d$igraph.obj, vertex.shape="rectangle", vertex.size=30 )
qgraph::qgraph( mod6d$adj.matrix )

#*****************************************************
# Model 7: Latent class analysis with 3 classes
#*****************************************************

#----  M7a: randomLCA (in randomLCA)
  #        - use two trials of starting values
mod7a <- randomLCA::randomLCA(dat,nclass=3, notrials=2, verbose=TRUE)
summary(mod7a)
plot(mod7a,type="l", xlab="Item")

#----  M7b: rasch.mirtlc (in sirt)
mod7b <- sirt::rasch.mirtlc( dat, Nclasses=3,seed=-30,  nstarts=2 )
summary(mod7b)
matplot( t(mod7b$pjk), type="l", xlab="Item" )

#----  M7c: poLCA (in poLCA)
  #   define formula for outcomes
f7c <- paste0( "cbind(", paste0(colnames(dat),collapse=","), ") ~ 1 " )
dat1 <- as.data.frame( dat + 1 ) # poLCA needs integer values from 1,2,..
mod7c <- poLCA::poLCA( stats::as.formula(f7c),dat1,nclass=3, verbose=TRUE)
plot(mod7c)

#----  M7d: gom.em (in sirt)
  #    - the latent class model is a special grade of membership model
mod7d <- sirt::gom.em( dat, K=3, problevels=c(0,1),  model="GOM"  )
summary(mod7d)

#---- - M7e: mirt (in mirt)
  # define three latent classes
Theta <- diag(3)
  # define mirt model
I <- ncol(dat)  # I=12
mirtmodel <- mirt::mirt.model("
        C1=1-12
        C2=1-12
        C3=1-12
        ")
  # get initial parameter values
mod.pars <- mirt::mirt(dat, model=mirtmodel,  pars="values")
  # modify parameters: only slopes refer to item-class probabilities
set.seed(9976)
  # set starting values for class specific item probabilities
mod.pars[ mod.pars$name=="d","value" ]  <- 0
mod.pars[ mod.pars$name=="d","est" ]  <- FALSE
b1 <- stats::qnorm( colMeans( dat ) )
mod.pars[ mod.pars$name=="a1","value" ]  <- b1
  # random starting values for other classes
mod.pars[ mod.pars$name \%in\% c("a2","a3"),"value" ]  <- b1 + stats::runif(12*2,-1,1)
mod.pars
  #** define prior for latent class analysis
lca_prior <- function(Theta,Etable){
  # number of latent Theta classes
  TP <- nrow(Theta)
  # prior in initial iteration
  if ( is.null(Etable) ){
    prior <- rep( 1/TP, TP )
  }
  # process Etable (this is correct for datasets without missing data)
  if ( ! is.null(Etable) ){
    # sum over correct and incorrect expected responses
    prior <- ( rowSums(Etable[, seq(1,2*I,2)]) + rowSums(Etable[,seq(2,2*I,2)]) )/I
  }
  prior <- prior / sum(prior)
  return(prior)
}
  #** estimate model
mod7e <- mirt::mirt(dat, mirtmodel, pars=mod.pars, verbose=TRUE,
            technical=list( customTheta=Theta, customPriorFun=lca_prior) )
  # compare estimated results
print(mod7e)
summary(mod7b)
  # The number of estimated parameters is incorrect because mirt does not correctly count
  # estimated parameters from the user customized  prior distribution.
mod7e@nest <- as.integer(sum(mod.pars$est) + 2)  # two additional class probabilities
  # extract log-likelihood
mod7e@logLik
  # compute AIC and BIC
( AIC <- -2*mod7e@logLik+2*mod7e@nest )
( BIC <- -2*mod7e@logLik+log(mod7e@Data$N)*mod7e@nest )
  # RMSEA and SRMSR fit statistic
mirt::M2(mod7e)     # TLI and CFI does not make sense in this example
  #** extract item parameters
sirt::mirt.wrapper.coef(mod7e)
  #** extract class-specific item-probabilities
probs <- apply( coef1[, c("a1","a2","a3") ], 2, stats::plogis )
matplot( probs, type="l", xlab="Item", main="mirt::mirt")
  #** inspect estimated distribution
mod7e@Theta
mod7e@Prior[[1]]

#*****************************************************
# Model 8: Mixed Rasch model with two classes
#*****************************************************

#----  M8a: raschmix (in psychomix)
mod8a <- psychomix::raschmix(data=as.matrix(dat), k=2, scores="saturated")
summary(mod8a)

#----  M8b: mrm (in mRm)
mod8b <- mRm::mrm(data.matrix=dat, cl=2)
mod8b$conv.to.bound
plot(mod8b)
print(mod8b)

#----  M8c: mirt (in mirt)
  #* define theta grid
theta.k <- seq( -5, 5, len=9 )
TP <- length(theta.k)
Theta <- matrix( 0, nrow=2*TP, ncol=4)
Theta[1:TP,1:2] <- cbind(theta.k, 1 )
Theta[1:TP + TP,3:4] <- cbind(theta.k, 1 )
Theta
  # define model
I <- ncol(dat)  # I=12
mirtmodel <- mirt::mirt.model("
        F1a=1-12  # slope Class 1
        F1b=1-12  # difficulty Class 1
        F2a=1-12  # slope Class 2
        F2b=1-12  # difficulty Class 2
        CONSTRAIN=(1-12,a1),(1-12,a3)
        ")
  # get initial parameter values
mod.pars <- mirt::mirt(dat, model=mirtmodel,  pars="values")
  # set starting values for class specific item probabilities
mod.pars[ mod.pars$name=="d","value" ]  <- 0
mod.pars[ mod.pars$name=="d","est" ]  <- FALSE
mod.pars[ mod.pars$name=="a1","value" ]  <- 1
mod.pars[ mod.pars$name=="a3","value" ]  <- 1
  # initial values difficulties
b1 <-  stats::qlogis( colMeans(dat) )
mod.pars[ mod.pars$name=="a2","value" ]  <- b1
mod.pars[ mod.pars$name=="a4","value" ]  <- b1 + stats::runif(I, -1, 1)
  #* define prior for mixed Rasch analysis
mixed_prior <- function(Theta,Etable){
  NC <- 2   # number of theta classes
  TP <- nrow(Theta) / NC
  prior1 <- stats::dnorm( Theta[1:TP,1] )
  prior1 <- prior1 / sum(prior1)
  if ( is.null(Etable) ){   prior <- c( prior1, prior1 ) }
  if ( ! is.null(Etable) ){
    prior <- ( rowSums( Etable[, seq(1,2*I,2)] ) +
                   rowSums( Etable[,seq(2,2*I,2)]) )/I
    a1 <- stats::aggregate( prior, list( rep(1:NC, each=TP) ), sum )
    a1[,2] <- a1[,2] / sum( a1[,2])
    # print some information during estimation
    cat( paste0( " Class proportions: ",
              paste0( round(a1[,2], 3 ), collapse=" " ) ), "\n")
    a1 <- rep( a1[,2], each=TP )
    # specify mixture of two normal distributions
    prior <- a1*c(prior1,prior1)
  }
  prior <- prior / sum(prior)
  return(prior)
}
  #* estimate model
mod8c <- mirt::mirt(dat, mirtmodel, pars=mod.pars, verbose=TRUE,
        technical=list(  customTheta=Theta, customPriorFun=mixed_prior ) )
  # Like in Model 7e, the number of estimated parameters must be included.
mod8c@nest <- as.integer(sum(mod.pars$est) + 1)
      # two class proportions and therefore one probability is freely estimated.
  #* extract item parameters
sirt::mirt.wrapper.coef(mod8c)
  #* estimated distribution
mod8c@Theta
mod8c@Prior

#----  M8d: tamaan (in TAM)

tammodel <- "
ANALYSIS:
  TYPE=MIXTURE ;
  NCLASSES(2);
  NSTARTS(7,20);
LAVAAN MODEL:
  F=~ A1__C4
  F ~~ F
ITEM TYPE:
  ALL(Rasch);
    "
mod8d <- TAM::tamaan( tammodel, resp=dat )
summary(mod8d)
# plot item parameters
I <- 12
ipars <- mod8d$itempartable_MIXTURE[ 1:I, ]
plot( 1:I, ipars[,3], type="o", ylim=range( ipars[,3:4] ), pch=16,
        xlab="Item", ylab="Item difficulty")
lines( 1:I, ipars[,4], type="l", col=2, lty=2)
points( 1:I, ipars[,4],  col=2, pch=2)

#*****************************************************
# Model 9: Mixed 2PL model with two classes
#*****************************************************

#----  M9a: tamaan (in TAM)

tammodel <- "
ANALYSIS:
  TYPE=MIXTURE ;
  NCLASSES(2);
  NSTARTS(10,30);
LAVAAN MODEL:
  F=~ A1__C4
  F ~~ F
ITEM TYPE:
  ALL(2PL);
    "
mod9a <- TAM::tamaan( tammodel, resp=dat )
summary(mod9a)

#*****************************************************
# Model 10: Rasch testlet model
#*****************************************************

#----  M10a: tam.fa (in TAM)
dims <- substring( colnames(dat),1,1 )  # define dimensions
mod10a <- TAM::tam.fa( resp=dat, irtmodel="bifactor1", dims=dims,
                control=list(maxiter=60) )
summary(mod10a)

#----  M10b: mirt (in mirt)
cmodel <- mirt::mirt.model("
        G=1-12
        A=1-4
        B=5-8
        C=9-12
        CONSTRAIN=(1-12,a1), (1-4, a2), (5-8, a3), (9-12,a4)
      ")
mod10b <- mirt::mirt(dat, model=cmodel, verbose=TRUE)
summary(mod10b)
coef(mod10b)
mod10b@logLik   # equivalent is slot( mod10b, "logLik")

#alternatively, using a dimensional reduction approach (faster and better accuracy)
cmodel <- mirt::mirt.model("
      G=1-12
      CONSTRAIN=(1-12,a1), (1-4, a2), (5-8, a3), (9-12,a4)
     ")
item_bundles <- rep(c(1,2,3), each=4)
mod10b1 <- mirt::bfactor(dat, model=item_bundles, model2=cmodel, verbose=TRUE)
coef(mod10b1)

#----  M10c: smirt (in sirt)
  # define Q-matrix
Qmatrix <- matrix(0,12,4)
Qmatrix[,1] <- 1
Qmatrix[ cbind( 1:12, match( dims, unique(dims)) +1 ) ]  <- 1
  # uncorrelated factors
variance.fixed <- cbind( c(1,1,1,2,2,3), c(2,3,4,3,4,4), 0 )
  # estimate model
mod10c <- sirt::smirt( dat, Qmatrix=Qmatrix, irtmodel="comp",
              variance.fixed=variance.fixed, qmcnodes=1000, maxiter=60)
summary(mod10c)

#*****************************************************
# Model 11: Bifactor model
#*****************************************************

#----  M11a: tam.fa (in TAM)
dims <- substring( colnames(dat),1,1 )  # define dimensions
mod11a <- TAM::tam.fa( resp=dat, irtmodel="bifactor2", dims=dims,
                 control=list(maxiter=60) )
summary(mod11a)

#----  M11b: bfactor (in mirt)
dims1 <- match( dims, unique(dims) )
mod11b <- mirt::bfactor(dat, model=dims1, verbose=TRUE)
summary(mod11b)
coef(mod11b)
mod11b@logLik

#----  M11c: smirt (in sirt)
  # define Q-matrix
Qmatrix <- matrix(0,12,4)
Qmatrix[,1] <- 1
Qmatrix[ cbind( 1:12, match( dims, unique(dims)) +1 ) ]  <- 1
  # uncorrelated factors
variance.fixed <- cbind( c(1,1,1,2,2,3), c(2,3,4,3,4,4), 0 )
  # estimate model
mod11c <- sirt::smirt( dat, Qmatrix=Qmatrix, irtmodel="comp", est.a="2PL",
                variance.fixed=variance.fixed, qmcnodes=1000, maxiter=60)
summary(mod11c)

#*****************************************************
# Model 12: Located latent class model: Rasch model with three theta classes
#*****************************************************

# use 10th item as the reference item
ref.item <- 10
# ability grid
theta.k <- seq(-4,4,len=9)

#----  M12a: rasch.mirtlc (in sirt)
mod12a <- sirt::rasch.mirtlc(dat, Nclasses=3, modeltype="MLC1", ref.item=ref.item)
summary(mod12a)

#----  M12b: gdm (in CDM)
theta.k <- seq(-1, 1, len=3)      # initial matrix
b.constraint <- matrix( c(10,1,0), nrow=1,ncol=3)
  # estimate model
mod12b <- CDM::gdm( dat, theta.k=theta.k, skillspace="est", irtmodel="1PL",
              b.constraint=b.constraint, maxiter=200)
summary(mod12b)

#----  M12c: mirt (in mirt)
items <- colnames(dat)
  # define three latent classes
Theta <- diag(3)
  # define mirt model
I <- ncol(dat)  # I=12
mirtmodel <- mirt::mirt.model("
        C1=1-12
        C2=1-12
        C3=1-12
        CONSTRAIN=(1-12,a1),(1-12,a2),(1-12,a3)
        ")
  # get parameters
mod.pars <- mirt(dat, model=mirtmodel,  pars="values")
 # set starting values for class specific item probabilities
mod.pars[ mod.pars$name=="d","value" ]  <- stats::qlogis( colMeans(dat,na.rm=TRUE) )
  # set item difficulty of reference item to zero
ind <- which( ( paste(mod.pars$item)==items[ref.item] ) &
               ( ( paste(mod.pars$name)=="d" ) ) )
mod.pars[ ind,"value" ]  <- 0
mod.pars[ ind,"est" ]  <- FALSE
  # initial values for a1, a2 and a3
mod.pars[ mod.pars$name \%in\% c("a1","a2","a3"),"value" ]  <- c(-1,0,1)
mod.pars
  #* define prior for latent class analysis
lca_prior <- function(Theta,Etable){
  # number of latent Theta classes
  TP <- nrow(Theta)
  # prior in initial iteration
  if ( is.null(Etable) ){
    prior <- rep( 1/TP, TP )
              }
  # process Etable (this is correct for datasets without missing data)
  if ( ! is.null(Etable) ){
    # sum over correct and incorrect expected responses
    prior <- ( rowSums( Etable[, seq(1,2*I,2)] ) + rowSums( Etable[, seq(2,2*I,2)] ) )/I
            }
  prior <- prior / sum(prior)
  return(prior)
   }
 #* estimate model
mod12c <- mirt(dat, mirtmodel, technical=list(
            customTheta=Theta, customPriorFun=lca_prior),
            pars=mod.pars, verbose=TRUE )
  # estimated parameters from the user customized  prior distribution.
mod12c@nest <- as.integer(sum(mod.pars$est) + 2)
  #* extract item parameters
coef1 <- sirt::mirt.wrapper.coef(mod12c)
  #* inspect estimated distribution
mod12c@Theta
coef1$coef[1,c("a1","a2","a3")]
mod12c@Prior[[1]]

#*****************************************************
# Model 13: Multidimensional model with discrete traits
#*****************************************************
# define Q-Matrix
Q <- matrix( 0, nrow=12,ncol=3)
Q[1:4,1] <- 1
Q[5:8,2] <- 1
Q[9:12,3] <- 1
# define discrete theta distribution with 3 dimensions
Theta <- scan(what="character",nlines=1)
  000 100 010 001 110 101 011 111
Theta <- as.numeric( unlist( lapply( Theta, strsplit, split="")   ) )
Theta <- matrix(Theta, 8, 3, byrow=TRUE )
Theta

#----  Model 13a: din (in CDM)
mod13a <- CDM::din( dat, q.matrix=Q, rule="DINA")
summary(mod13a)
# compare used Theta distributions
cbind( Theta, mod13a$attribute.patt.splitted)

#----  Model 13b: gdm (in CDM)
mod13b <- CDM::gdm( dat, Qmatrix=Q, theta.k=Theta, skillspace="full")
summary(mod13b)

#----  Model 13c: mirt (in mirt)
  # define mirt model
I <- ncol(dat)  # I=12
mirtmodel <- mirt::mirt.model("
        F1=1-4
        F2=5-8
        F3=9-12
        ")
  # get parameters
mod.pars <- mirt(dat, model=mirtmodel,  pars="values")
# starting values d parameters (transformed guessing parameters)
ind <- which(  mod.pars$name=="d"  )
mod.pars[ind,"value"] <- stats::qlogis(.2)
# starting values transformed slipping parameters
ind <- which( ( mod.pars$name \%in\% paste0("a",1:3)  ) &  ( mod.pars$est ) )
mod.pars[ind,"value"] <- stats::qlogis(.8) - stats::qlogis(.2)
mod.pars

  #* define prior for latent class analysis
lca_prior <- function(Theta,Etable){
  TP <- nrow(Theta)
  if ( is.null(Etable) ){
    prior <- rep( 1/TP, TP )
              }
  if ( ! is.null(Etable) ){
    prior <- ( rowSums( Etable[, seq(1,2*I,2)] ) + rowSums( Etable[, seq(2,2*I,2)] ) )/I
            }
  prior <- prior / sum(prior)
  return(prior)
}
 #* estimate model
mod13c <- mirt(dat, mirtmodel, technical=list(
            customTheta=Theta, customPriorFun=lca_prior),
            pars=mod.pars, verbose=TRUE )
  # estimated parameters from the user customized  prior distribution.
mod13c@nest <- as.integer(sum(mod.pars$est) + 2)
  #* extract item parameters
coef13c <- sirt::mirt.wrapper.coef(mod13c)$coef
  #* inspect estimated distribution
mod13c@Theta
mod13c@Prior[[1]]

 #-* comparisons of estimated  parameters
# extract guessing and slipping parameters from din
dfr <- coef(mod13a)[, c("guess","slip") ]
colnames(dfr) <- paste0("din.",c("guess","slip") )
# estimated parameters from gdm
dfr$gdm.guess <- stats::plogis(mod13b$item$b)
dfr$gdm.slip <- 1 - stats::plogis( rowSums(mod13b$item[,c("b.Cat1","a.F1","a.F2","a.F3")] ) )
# estimated parameters from mirt
dfr$mirt.guess <- stats::plogis( coef13c$d )
dfr$mirt.slip <- 1 - stats::plogis( rowSums(coef13c[,c("d","a1","a2","a3")]) )
# comparison
round(dfr[, c(1,3,5,2,4,6)],3)
  ##      din.guess gdm.guess mirt.guess din.slip gdm.slip mirt.slip
  ##   A1     0.691     0.684      0.686    0.000    0.000     0.000
  ##   A2     0.491     0.489      0.489    0.031    0.038     0.036
  ##   A3     0.302     0.300      0.300    0.184    0.193     0.190
  ##   A4     0.244     0.239      0.240    0.337    0.340     0.339
  ##   B1     0.568     0.579      0.577    0.163    0.148     0.151
  ##   B2     0.329     0.344      0.340    0.344    0.326     0.329
  ##   B3     0.817     0.827      0.825    0.014    0.007     0.009
  ##   B4     0.431     0.463      0.456    0.104    0.089     0.092
  ##   C1     0.188     0.191      0.189    0.013    0.013     0.013
  ##   C2     0.050     0.050      0.050    0.239    0.238     0.239
  ##   C3     0.000     0.002      0.001    0.065    0.065     0.065
  ##   C4     0.000     0.004      0.000    0.212    0.212     0.212

# estimated class sizes
dfr <- data.frame( "Theta"=Theta, "din"=mod13a$attribute.patt$class.prob,
                   "gdm"=mod13b$pi.k, "mirt"=mod13c@Prior[[1]])
# comparison
round(dfr,3)
  ##     Theta.1 Theta.2 Theta.3   din   gdm  mirt
  ##   1       0       0       0 0.039 0.041 0.040
  ##   2       1       0       0 0.008 0.009 0.009
  ##   3       0       1       0 0.009 0.007 0.008
  ##   4       0       0       1 0.394 0.417 0.412
  ##   5       1       1       0 0.011 0.011 0.011
  ##   6       1       0       1 0.017 0.042 0.037
  ##   7       0       1       1 0.042 0.008 0.016
  ##   8       1       1       1 0.480 0.465 0.467

#*****************************************************
# Model 14: DINA model with two skills
#*****************************************************

# define some simple Q-Matrix (does not really make in this application)
Q <- matrix( 0, nrow=12,ncol=2)
Q[1:4,1] <- 1
Q[5:8,2] <- 1
Q[9:12,1:2] <- 1
# define discrete theta distribution with 3 dimensions
Theta <- scan(what="character",nlines=1)
  00 10 01 11
Theta <- as.numeric( unlist( lapply( Theta, strsplit, split="")   ) )
Theta <- matrix(Theta, 4, 2, byrow=TRUE )
Theta

#----  Model 14a: din (in CDM)
mod14a <- CDM::din( dat, q.matrix=Q, rule="DINA")
summary(mod14a)
# compare used Theta distributions
cbind( Theta, mod14a$attribute.patt.splitted)

#----  Model 14b: mirt (in mirt)
  # define mirt model
I <- ncol(dat)  # I=12
mirtmodel <- mirt::mirt.model("
        F1=1-4
        F2=5-8
        (F1*F2)=9-12
        ")
#-> constructions like (F1*F2*F3) are also allowed in mirt.model
  # get parameters
mod.pars <- mirt(dat, model=mirtmodel,  pars="values")
# starting values d parameters (transformed guessing parameters)
ind <- which(  mod.pars$name=="d"  )
mod.pars[ind,"value"] <- stats::qlogis(.2)
# starting values transformed slipping parameters
ind <- which( ( mod.pars$name \%in\% paste0("a",1:3)  ) &  ( mod.pars$est ) )
mod.pars[ind,"value"] <- stats::qlogis(.8) - stats::qlogis(.2)
mod.pars
 #* use above defined prior lca_prior
 # lca_prior <- function(prior,Etable) ...
 #* estimate model
mod14b <- mirt(dat, mirtmodel, technical=list(
            customTheta=Theta, customPriorFun=lca_prior),
            pars=mod.pars, verbose=TRUE )
  # estimated parameters from the user customized  prior distribution.
mod14b@nest <- as.integer(sum(mod.pars$est) + 2)
  #* extract item parameters
coef14b <- sirt::mirt.wrapper.coef(mod14b)$coef

 #-* comparisons of estimated  parameters
# extract guessing and slipping parameters from din
dfr <- coef(mod14a)[, c("guess","slip") ]
colnames(dfr) <- paste0("din.",c("guess","slip") )
# estimated parameters from mirt
dfr$mirt.guess <- stats::plogis( coef14b$d )
dfr$mirt.slip <- 1 - stats::plogis( rowSums(coef14b[,c("d","a1","a2","a3")]) )
# comparison
round(dfr[, c(1,3,2,4)],3)
  ##      din.guess mirt.guess din.slip mirt.slip
  ##   A1     0.674      0.671    0.030     0.030
  ##   A2     0.423      0.420    0.049     0.050
  ##   A3     0.258      0.255    0.224     0.225
  ##   A4     0.245      0.243    0.394     0.395
  ##   B1     0.534      0.543    0.166     0.164
  ##   B2     0.338      0.347    0.382     0.380
  ##   B3     0.796      0.802    0.016     0.015
  ##   B4     0.421      0.436    0.142     0.140
  ##   C1     0.850      0.851    0.000     0.000
  ##   C2     0.480      0.480    0.097     0.097
  ##   C3     0.746      0.746    0.026     0.026
  ##   C4     0.575      0.577    0.136     0.137

# estimated class sizes
dfr <- data.frame( "Theta"=Theta, "din"=mod13a$attribute.patt$class.prob,
                    "mirt"=mod14b@Prior[[1]])
# comparison
round(dfr,3)
  ##     Theta.1 Theta.2   din  mirt
  ##   1       0       0 0.357 0.369
  ##   2       1       0 0.044 0.049
  ##   3       0       1 0.047 0.031
  ##   4       1       1 0.553 0.551

#*****************************************************
# Model 15: Rasch model with non-normal distribution
#*****************************************************

# A non-normal theta distributed is specified by log-linear smoothing
# the distribution as described in
# Xu, X., & von Davier, M. (2008). Fitting the structured general diagnostic model
# to NAEP data. ETS Research Report ETS RR-08-27. Princeton, ETS.

# define theta grid
theta.k <- matrix( seq(-4,4,len=15), ncol=1 )
# define design matrix for smoothing (up to cubic moments)
delta.designmatrix <- cbind( 1, theta.k, theta.k^2, theta.k^3 )
# constrain item difficulty of fifth item (item B1) to zero
b.constraint <- matrix( c(5,1,0), ncol=3 )

#----  Model 15a: gdm (in CDM)
mod15a <- CDM::gdm( dat, irtmodel="1PL", theta.k=theta.k,
               b.constraint=b.constraint  )
summary(mod15a)
 # plot estimated distribution
graphics::barplot( mod15a$pi.k[,1], space=0, names.arg=round(theta.k[,1],2),
           main="Estimated Skewed Distribution (gdm function)")

#----  Model 15b: mirt (in mirt)
 # define mirt model
mirtmodel <- mirt::mirt.model("
    F=1-12
    ")
 # get parameters
mod.pars <- mirt::mirt(dat, model=mirtmodel, pars="values", itemtype="Rasch")
  # fix variance (just for correct counting of parameters)
mod.pars[ mod.pars$name=="COV_11", "est"] <- FALSE
  # fix item difficulty
ind <- which( ( mod.pars$item=="B1" ) & ( mod.pars$name=="d" ) )
mod.pars[ ind, "value"] <- 0
mod.pars[ ind, "est"] <- FALSE

 # define prior
loglinear_prior <- function(Theta,Etable){
    TP <- nrow(Theta)
    if ( is.null(Etable) ){
    prior <- rep( 1/TP, TP )
           }
    # process Etable (this is correct for datasets without missing data)
    if ( ! is.null(Etable) ){
          # sum over correct and incorrect expected responses
       prior <- ( rowSums( Etable[, seq(1,2*I,2)] ) + rowSums( Etable[, seq(2,2*I,2)] ) )/I
       # smooth prior using the above design matrix and a log-linear model
       # see Xu & von Davier (2008).
       y <- log( prior + 1E-15 )
       lm1 <- lm( y ~ 0 + delta.designmatrix, weights=prior )
       prior <- exp(fitted(lm1))   # smoothed prior
           }
    prior <- prior / sum(prior)
    return(prior)
}

#* estimate model
mod15b <- mirt::mirt(dat, mirtmodel, technical=list(
                customTheta=theta.k, customPriorFun=loglinear_prior ),
                pars=mod.pars, verbose=TRUE )
# estimated parameters from the user customized prior distribution.
mod15b@nest <- as.integer(sum(mod.pars$est) + 3)
#* extract item parameters
coef1 <- sirt::mirt.wrapper.coef(mod15b)$coef

#** compare estimated item parameters
dfr <- data.frame( "gdm"=mod15a$item$b.Cat1, "mirt"=coef1$d )
rownames(dfr) <- colnames(dat)
round(t(dfr),4)
  ##            A1     A2      A3      A4 B1      B2     B3     B4     C1    C2     C3    C4
  ##   gdm  0.9818 0.1538 -0.7837 -1.3197  0 -1.0902 1.6088 -0.170 1.9778 0.006 1.1859 0.135
  ##   mirt 0.9829 0.1548 -0.7826 -1.3186  0 -1.0892 1.6099 -0.169 1.9790 0.007 1.1870 0.136
# compare estimated theta distribution
dfr <- data.frame( "gdm"=mod15a$pi.k, "mirt"=mod15b@Prior[[1]] )
round(t(dfr),4)
  ##        1 2     3     4      5      6      7      8      9     10     11     12     13
  ##   gdm  0 0 1e-04 9e-04 0.0056 0.0231 0.0652 0.1299 0.1881 0.2038 0.1702 0.1129 0.0612
  ##   mirt 0 0 1e-04 9e-04 0.0056 0.0232 0.0653 0.1300 0.1881 0.2038 0.1702 0.1128 0.0611
  ##            14    15
  ##   gdm  0.0279 0.011
  ##   mirt 0.0278 0.011
}
}

\keyword{datasets}
