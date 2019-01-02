%% File Name: R2conquest.Rd
%% File Version: 3.23

\name{R2conquest}
\alias{R2conquest}
\alias{summary.R2conquest}
\alias{read.show}
\alias{read.show.term}
\alias{read.show.regression}
\alias{read.pv}
\alias{read.multidimpv}
\alias{read.pimap}



%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Running ConQuest From Within \R
}

\description{
The function \code{R2conquest} runs the IRT software ConQuest
(Wu, Adams, Wilson & Haldane, 2007) from within \R.

Other functions are utility functions for reading item parameters, plausible
values or person-item maps.
}
\usage{
R2conquest(dat, path.conquest, conquest.name="console", converge=0.001,
    deviancechange=1e-04, iter=800, nodes=20, minnode=-6, maxnode=6,
    show.conquestoutput=FALSE, name="rasch", pid=1:(nrow(dat)), wgt=NULL, X=NULL,
    set.constraints=NULL, model="item", regression=NULL,
    itemcodes=seq(0,max(dat,na.rm=TRUE)), constraints=NULL, digits=5, onlysyntax=FALSE,
    qmatrix=NULL, import.regression=NULL, anchor.regression=NULL,
    anchor.covariance=NULL, pv=TRUE, designmatrix=NULL, only.calibration=FALSE,
    init_parameters=NULL, n_plausible=10,  persons.elim=TRUE, est.wle=TRUE,
    save.bat=TRUE, use.bat=FALSE, read.output=TRUE, ignore.pid=FALSE)

\method{summary}{R2conquest}(object, \dots)

# read all terms in a show file or only some terms
read.show(showfile)
read.show.term(showfile, term)

# read regression parameters in a show file
read.show.regression(showfile)

# read unidimensional plausible values form a pv file
read.pv(pvfile, npv=5)
# read multidimensional plausible values
read.multidimpv(pvfile, ndim, npv=5)

# read person-item map
read.pimap(showfile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame of item responses
}
  \item{path.conquest}{
Directory where the ConQuest executable file is located
}
  \item{conquest.name}{
Name of the ConQuest executable.
}
  \item{converge}{
Maximal change in parameters
}
  \item{deviancechange}{
Maximal change in deviance
}
  \item{iter}{
Maximum number of iterations
}
  \item{nodes}{
Number of nodes for integration
}
\item{minnode}{Minimum value of discrete grid of \eqn{\theta} nodes}
\item{maxnode}{Maximum value of discrete grid of \eqn{\theta} nodes}
  \item{show.conquestoutput}{
Show ConQuest run log file on console?
}
  \item{name}{
Name of the output files. The default is \code{'rasch'}.
}
  \item{pid}{
Person identifier
}
  \item{wgt}{
Vector of person weights
}
  \item{X}{
Matrix of covariates for the latent regression model
(e.g. gender, socioeconomic status, ..) or
for the item design (e.g. raters, booklets, ...)
}
  \item{set.constraints}{
This is the set.constraints in ConQuest. It can be
\code{"cases"} (constraint for persons), \code{"items"}
or \code{"none"}
}
  \item{model}{
Definition model statement. It can be for example
\code{"item+item*step"} or \code{"item+booklet+rater"}
}
  \item{regression}{
The ConQuest regression statement (for example \code{"gender+status"})
}
  \item{itemcodes}{
Vector of valid codes for item responses. E.g. for partial credit
data with at most 3 points it must be \code{c(0,1,2,3)}.
}
  \item{constraints}{
Matrix of item parameter constraints. 1st column: Item names,
2nd column: Item parameters. It only works correctly
for dichotomous data.
}
  \item{digits}{
Number of digits for covariates in the latent regression model
}
  \item{onlysyntax}{
Should only be ConQuest syntax generated?
}
  \item{qmatrix}{
Matrix of item loadings on dimensions in a multidimensional
IRT model
}
  \item{import.regression}{
Name of an file with initial covariance parameters
(follow the ConQuest specification rules!)
}
  \item{anchor.regression}{
Name of an file with anchored regression parameters
}
  \item{anchor.covariance}{
Name of an file with anchored covariance parameters
(follow the ConQuest specification rules!)
}
  \item{pv}{
Draw plausible values?
}
  \item{designmatrix}{
Design matrix for item parameters (see the ConQuest manual)
}
  \item{only.calibration}{
Estimate only item parameters and not person parameters
(no WLEs or plausible values are estimated)?
}
  \item{init_parameters}{
Name of an file with initial item parameters (follow the
ConQuest specification rules!)
}
  \item{n_plausible}{
Number of plausible values
}
  \item{persons.elim}{
Eliminate persons with only missing item responses?
}
  \item{est.wle}{
Estimate weighted likelihood estimate?
}
  \item{save.bat}{
Save bat file?
}
  \item{use.bat}{
Run ConQuest from within \R due a direct call via the \code{system}
command (\code{use.bat=FALSE}) or via a \code{system} call
of a bat file in the working directory (\code{use.bat=TRUE})
}
  \item{read.output}{
Should ConQuest output files be processed? Default is \code{TRUE}.
}
\item{ignore.pid}{Logical indicating whether person identifiers (\code{pid})
should be processed in ConQuest input syntax.}
\item{object}{
Object of class \code{R2conquest}
}
\item{showfile}{
A ConQuest show file (\code{shw} file)
}
\item{term}{
Name of the term to be extracted in the show file
}
\item{pvfile}{File with plausible values}
\item{ndim}{Number of dimensions}
\item{npv}{Number of plausible values}
\item{\dots}{Further arguments to be passed}
}


\details{
Consult the ConQuest manual (Wu et al., 2007) for specification
details.
}



\value{
A list with several entries

\item{item}{Data frame with item parameters and item statistics}
\item{person}{Data frame with person parameters}
\item{shw.itemparameter}{ConQuest output table for item parameters}
\item{shw.regrparameter}{ConQuest output table for regression parameters}
\item{\dots}{More values}
}


\references{
Wu, M. L., Adams, R. J., Wilson, M. R. & Haldane, S. (2007).
\emph{ACER ConQuest Version 2.0}. Mulgrave.
{https://shop.acer.edu.au/acer-shop/group/CON3}.
}


%\author{
%Alexander Robitzsch
%}



%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also the \pkg{eat} package (\url{http://r-forge.r-project.org/projects/eat/})
for elaborate functionality of using ConQuest from within \R.

See also the \pkg{TAM} package for similar (and even extended)
functionality for specifying item response models.
}
\examples{
\dontrun{
# define ConQuest path
path.conquest <- "C:/Conquest/"

#############################################################################
# EXAMPLE 1: Dichotomous data (data.pisaMath)
#############################################################################
library(sirt)
data(data.pisaMath)
dat <- data.pisaMath$data

# select items
items <- colnames(dat)[ which( substring( colnames(dat), 1, 1)=="M" ) ]

#***
# Model 11: Rasch model
mod11 <- sirt::R2conquest(dat=dat[,items], path.conquest=path.conquest,
             pid=dat$idstud, name="mod11")
summary(mod11)
# read show file
shw11 <- sirt::read.show( "mod11.shw" )
# read person-item map
pi11 <- sirt::read.pimap(showfile="mod11.shw")

#***
# Model 12: Rasch model with fixed item difficulties (from Model 1)
mod12 <- sirt::R2conquest(dat=dat[,items], path.conquest=path.conquest,
             pid=dat$idstud, constraints=mod11$item[, c("item","itemdiff")],
             name="mod12")
summary(mod12)

#***
# Model 13: Latent regression model with predictors female, hisei and migra
mod13a <- sirt::R2conquest(dat=dat[,items], path.conquest=path.conquest,
             pid=dat$idstud, X=dat[, c("female", "hisei", "migra") ],
             name="mod13a")
summary(mod13a)

# latent regression with a subset of predictors
mod13b <- sirt::R2conquest(dat=dat[,items], path.conquest=path.conquest,
             pid=dat$idstud, X=dat[, c("female", "hisei", "migra") ],
             regression="hisei migra", name="mod13b")

#***
# Model 14: Differential item functioning (female)
mod14 <- sirt::R2conquest(dat=dat[,items], path.conquest=path.conquest,
             pid=dat$idstud, X=dat[, c("female"), drop=FALSE],
             model="item+female+item*female",  regression="",  name="mod14")

#############################################################################
# EXAMPLE 2: Polytomous data (data.Students)
#############################################################################
library(CDM)
data(data.Students)
dat <- data.Students

# select items
items <- grep.vec( "act", colnames(dat) )$x

#***
# Model 21: Partial credit model
mod21 <- sirt::R2conquest(dat=dat[,items], path.conquest=path.conquest,
              model="item+item*step",  name="mod21")

#***
# Model 22: Rating scale model
mod22 <- sirt::R2conquest(dat=dat[,items], path.conquest=path.conquest,
              model="item+step", name="mod22")

#***
# Model 23: Multidimensional model
items <- grep.vec( c("act", "sc" ), colnames(dat),  "OR" )$x
qmatrix <- matrix( 0, nrow=length(items), 2 )
qmatrix[1:5,1] <- 1
qmatrix[6:9,2] <- 1
mod23 <- sirt::R2conquest(dat=dat[,items], path.conquest=path.conquest,
            model="item+item*step", qmatrix=qmatrix, name="mod23")

#############################################################################
# EXAMPLE 3: Multi facet models (data.ratings1)
#############################################################################
library(sirt)
data(data.ratings1)
dat <- data.ratings1

items <- paste0("k",1:5)

# use numeric rater ID's
raters <- as.numeric( substring( paste( dat$rater ), 3 ) )

#***
# Model 31: Rater model 'item+item*step+rater'
mod31 <- sirt::R2conquest(dat=dat[,items], path.conquest=path.conquest,
              itemcodes=0:3, model="item+item*step+rater",
              pid=dat$idstud, X=data.frame("rater"=raters),
              regression="", name="mod31")

#***
# Model 32: Rater model 'item+item*step+rater+item*rater'
mod32 <- sirt::R2conquest(dat=dat[,items], path.conquest=path.conquest,
              model="item+item*step+rater+item*rater",
              pid=dat$idstud, X=data.frame("rater"=raters),
              regression="", name="mod32")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ConQuest}
%% \keyword{summary}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
