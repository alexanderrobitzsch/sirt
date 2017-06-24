\name{data.eid}
\alias{data.eid}
\docType{data}
\title{
Examples with Datasets from Eid and Schmidt (2014)
}

\description{
Examples with datasets from Eid and Schmidt (2014),
illustrations with several \R packages. The examples
follow closely the online material of Hosoya (2014).
}

\usage{
data(data.eid)
}

\format{
The dataset \code{data.eid} is just a placeholder.
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
The material can be downloaded from
\url{http://www.hogrefe.de/buecher/lehrbuecher/psychlehrbuchplus/lehrbuecher/
testtheorie-und-testkonstruktion/zusatzmaterial/}.
}

\references{
Eid, M., & Schmidt, K. (2014). \emph{Testtheorie und Testkonstruktion}.
Goettingen, Hogrefe.

Hosoya, G. (2014). \emph{Einfuehrung in die Analyse testtheoretischer
Modelle mit \R}. Available at 
\emph{http://www.hogrefe.de/buecher/lehrbuecher/psychlehrbuchplus/lehrbuecher/testtheorie-und-testkonstruktion/zusatzmaterial/}.
}


\examples{
\dontrun{	
# The "dataset" data.eid is just a placeholder.
site <- paste0( "http://www.hogrefe.de/fileadmin/redakteure/hogrefe_de/" ,
        "Psychlehrbuchplus/Testtheorie_und_Testkonstruktion/R-Analysen/" )

miceadds::library_install("foreign")
#---- load some IRT packages in R
miceadds::library_install("TAM")        # package (a)
miceadds::library_install("mirt")       # package (b)
miceadds::library_install("sirt")       # package (c)
miceadds::library_install("eRm")        # package (d)
miceadds::library_install("ltm")        # package (e)
miceadds::library_install("psychomix")  # package (f)

#############################################################################
# EXAMPLES Ch. 4: Unidimensional IRT models | dichotomous data
#############################################################################

# link to dataset
linkname <- paste0( site , "ids_new.sav")
# load data
data0 <- foreign::read.spss( linkname , to.data.frame=TRUE , use.value.labels=FALSE)
# extract items
dat <- data0[,2:11]

#*********************************************************
# Model 1: Rasch model
#*********************************************************

#-----------
#-- 1a: estimation with TAM package

# estimation with tam.mml
mod1a <- TAM::tam.mml(dat)
summary(mod1a)

# person parameters in TAM
pp1a <- TAM::tam.wle(mod1a)

# plot item response functions
plot(mod1a,export=FALSE,ask=TRUE)

# Infit and outfit in TAM
itemf1a <- TAM::tam.fit(mod1a)
itemf1a

# model fit
modf1a <- TAM::tam.modelfit(mod1a)
summary(modf1a)

#-----------
#-- 1b: estimation with mirt package

# estimation with mirt
mod1b <- mirt::mirt( dat , 1 , itemtype="Rasch")
summary(mod1b)
print(mod1b)

# person parameters
pp1b <- mirt::fscores(mod1b , method="WLE")

# extract coefficients
sirt::mirt.wrapper.coef(mod1b)

# plot item response functions
plot(mod1b, type="trace" )
par(mfrow=c(1,1))

# item fit
itemf1b <- mirt::itemfit(mod1b)
itemf1b

# model fit
modf1b <- mirt::M2(mod1b)
modf1b

#-----------
#-- 1c: estimation with sirt package

# estimation with rasch.mml2
mod1c <- sirt::rasch.mml2(dat)
summary(mod1c)

# person parameters (EAP)
pp1c <- mod1c$person

# plot item response functions
plot(mod1c , ask=TRUE )

# model fit
modf1c <- sirt::modelfit.sirt(mod1c)
summary(modf1c)

#-----------
#-- 1d: estimation with eRm package

# estimation with RM
mod1d <- eRm::RM(dat)
summary(mod1d)

# estimation person parameters
pp1d <- eRm::person.parameter(mod1d)
summary(pp1d)

# plot item response functions
eRm::plotICC(mod1d)

# person-item map
eRm::plotPImap(mod1d)

# item fit
itemf1d <- eRm::itemfit(pp1d)

# person fit
persf1d <- eRm::personfit(pp1d)

#-----------
#-- 1e: estimation with ltm package

# estimation with rasch
mod1e <- ltm::rasch(dat)
summary(mod1e)

# estimation person parameters
pp1e <- ltm::factor.scores(mod1e)

# plot item response functions
plot(mod1e)

# item fit
itemf1e <- ltm::item.fit(mod1e)

# person fit
persf1e <- ltm::person.fit(mod1e)

# goodness of fit with Bootstrap
modf1e <- ltm::GoF.rasch(mod1e,B=20)    # use more bootstrap samples
modf1e

#*********************************************************
# Model 2: 2PL model
#*********************************************************

#-----------
#-- 2a: estimation with TAM package

# estimation
mod2a <- TAM::tam.mml.2pl(dat)
summary(mod2a)

# model fit
modf2a <- TAM::tam.modelfit(mod2a)
summary(modf2a)

# item response functions
plot(mod2a , export=FALSE , ask=TRUE)

# model comparison
anova(mod1a,mod2a)

#-----------
#-- 2b: estimation with mirt package

# estimation
mod2b <- mirt::mirt(dat,1,itemtype="2PL")
summary(mod2b)
print(mod2b)
sirt::mirt.wrapper.coef(mod2b)

# model fit
modf2b <- mirt::M2(mod2b)
modf2b

#-----------
#-- 2c: estimation with sirt package

I <- ncol(dat)
# estimation
mod2c <- sirt::rasch.mml2(dat,est.a=1:I)
summary(mod2c)

# model fit
modf2c <- sirt::modelfit.sirt(mod2c)
summary(modf2c)

#-----------
#-- 2e: estimation with ltm package

# estimation
mod2e <- ltm::ltm(dat ~ z1 )
summary(mod2e)

# item response functions
plot(mod2e)

#*********************************************************
# Model 3: Mixture Rasch model
#*********************************************************

#-----------
#-- 3a: estimation with TAM package

# avoid "_" in column names if the "__" operator is used in 
# the tamaan syntax
dat1 <- dat
colnames(dat1) <- gsub("_" , "" , colnames(dat1) )
# define tamaan model
tammodel <- "
ANALYSIS:
  TYPE=MIXTURE ;
  NCLASSES(2);
  NSTARTS(20,25);   # 20 random starts with 25 initial iterations each
LAVAAN MODEL:
  F =~ Freude1__Freude2
  F ~~ F
ITEM TYPE:
  ALL(Rasch);
    "    
mod3a <- TAM::tamaan( tammodel , resp=dat1 )
summary(mod3a)
# extract item parameters
ipars <- mod2$itempartable_MIXTURE[ 1:10 , ]
plot( 1:10 , ipars[,3] , type="o" , ylim= range( ipars[,3:4] ) , pch=16 ,
        xlab="Item" , ylab="Item difficulty")
lines( 1:10 , ipars[,4] , type="l", col=2 , lty=2)
points( 1:10 , ipars[,4] ,  col=2 , pch=2)

#-----------
#-- 3f: estimation with psychomix package

# estimation
mod3f <- psychomix::raschmix( as.matrix(dat) , k=2 , scores="meanvar")
summary(mod3f)
# plot class-specific item difficulties
plot(mod3f)

#############################################################################
# EXAMPLES Ch. 5: Unidimensional IRT models | polytomous data
#############################################################################

# link to dataset
linkname <- paste0( site , "Daten-kapitel-5-sex.sav")
# load data
data0 <- foreign::read.spss( linkname , to.data.frame=TRUE , use.value.labels=FALSE)
# extract items
dat <- data0[,2:7]

#*********************************************************
# Model 1: Partial credit model
#*********************************************************

#-----------
#-- 1a: estimation with TAM package

# estimation with tam.mml
mod1a <- TAM::tam.mml(dat)
summary(mod1a)

# person parameters in TAM
pp1a <- tam.wle(mod1a)

# plot item response functions
plot(mod1a,export=FALSE,ask=TRUE)

# Infit and outfit in TAM
itemf1a <- TAM::tam.fit(mod1a)
itemf1a

# model fit
modf1a <- TAM::tam.modelfit(mod1a)
summary(modf1a)

#-----------
#-- 1b: estimation with mirt package

# estimation with tam.mml
mod1b <- mirt::mirt( dat , 1 , itemtype="Rasch")
summary(mod1b)
print(mod1b)
sirt::mirt.wrapper.coef(mod1b)

# plot item response functions
plot(mod1b, type="trace" )
par(mfrow=c(1,1))

# item fit
itemf1b <- mirt::itemfit(mod1b)
itemf1b

#-----------
#-- 1c: estimation with sirt package

# estimation with rm.facets
mod1c <- sirt::rm.facets(dat)
summary(mod1c)
summary(mod1a)

#-----------
#-- 1d: estimation with eRm package

# estimation
mod1d <- eRm::PCM(dat)
summary(mod1d)

# plot item response functions
eRm::plotICC(mod1d)

# person-item map
eRm::plotPImap(mod1d)

# item fit
itemf1d <- eRm::itemfit(pp1d)

#-----------
#-- 1e: estimation with ltm package

# estimation
mod1e <- ltm::gpcm(dat, constraint="1PL")
summary(mod1e)
# plot item response functions
plot(mod1e)

#*********************************************************
# Model 2: Generalized partial credit model
#*********************************************************

#-----------
#-- 2a: estimation with TAM package

# estimation with tam.mml
mod2a <- TAM::tam.mml.2pl(dat, irtmodel="GPCM")
summary(mod2a)

# model fit
modf2a <- TAM::tam.modelfit(mod2a)
summary(modf2a)

#-----------
#-- 2b: estimation with mirt package

# estimation
mod2b <- mirt::mirt( dat , 1 , itemtype="gpcm")
summary(mod2b)
print(mod2b)
sirt::mirt.wrapper.coef(mod2b)

#-----------
#-- 2c: estimation with sirt package

# estimation with rm.facets
mod2c <- sirt::rm.facets(dat , est.a.item=TRUE)
summary(mod2c)

#-----------
#-- 2e: estimation with ltm package

# estimation
mod2e <- ltm::gpcm(dat)
summary(mod2e)
plot(mod2e)
}
}


\keyword{datasets}
