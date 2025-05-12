%% File Name: xxirt_createParTable.Rd
%% File Version: 1.172

\name{xxirt_createParTable}
\alias{xxirt_createParTable}
\alias{xxirt_createDiscItem}
\alias{xxirt_modifyParTable}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create Item Response Functions and Item Parameter Table
}
\description{
Create item response functions and item parameter table
}

\usage{
xxirt_createDiscItem( name, par, est, P, lower=-Inf, upper=Inf,
     prior=NULL, prior_par1=NULL, prior_par2=NULL, X=NULL)

xxirt_createParTable(dat, itemtype, customItems=NULL)

xxirt_modifyParTable( partable, parname, item=NULL, value=NULL,
     est=NULL, parlabel=NULL, parindex=NULL, lower=NULL,
     upper=NULL, prior=NULL, prior_par1=NULL, prior_par2=NULL )
}



%- maybe also 'usage' for other objects documented here.
\arguments{
\item{name}{Type of item response function}
\item{par}{Named vector of starting values of item parameters}
\item{est}{Logical vector indicating which parameters should be estimated}
\item{P}{Item response function}
\item{lower}{Lower bounds}
\item{upper}{Upper bounds}
\item{prior}{Prior distribution}
\item{prior_par1}{First parameter prior distribution}
\item{prior_par2}{Second parameter prior distribution}
\item{X}{optional matrix of person covariates}
  \item{dat}{
Data frame with item responses
}
  \item{itemtype}{
Vector of item types
}
 \item{customItems}{
List with item objects created by \code{xxirt_createDiscItem}
}
\item{partable}{Item parameter table}
\item{parname}{Parameter name}
\item{item}{Item}
\item{value}{Value of item parameter}
\item{parindex}{Parameter index}
\item{parlabel}{Item parameter label}
}

% xxirt_modifyParTable( partable, parname, item=NULL, value=NULL,
%     est=NULL, parlabel=NULL, parindex=NULL, lower=NULL,
%     upper=NULL, prior=NULL, prior_par1=NULL, prior_par2=NULL )


%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}

%\author{
%Alexander Robitzsch
%}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{xxirt}}

See \code{\link[mirt:createItem]{mirt::createItem}} for similar
functionality.
}

\examples{
#############################################################################
## EXAMPLE 1: Definition of item response functions
#############################################################################

data(data.read)
dat <- data.read

#------ Definition of item response functions
#*** IRF 2PL
P_2PL <- function( par, Theta, ncat){
    a <- par[1]
    b <- par[2]
    TP <- nrow(Theta)
    P <- matrix( NA, nrow=TP, ncol=ncat)
    P[,1] <- 1
    for (cc in 2:ncat){
        P[,cc] <- exp( (cc-1) * a * Theta[,1] - b )
    }
    P <- P / rowSums(P)
    return(P)
}

#*** IRF 1PL
P_1PL <- function( par, Theta, ncat){
    b <- par[1]
    TP <- nrow(Theta)
    par0 <- c(1,b)
    P <- P_2PL( par=par0, Theta=Theta, ncat=ncat)
    return(P)
}

#** created item classes of 1PL and 2PL models
par <- c( "a"=1, "b"=0 )
# define some slightly informative prior of 2PL
item_2PL <- sirt::xxirt_createDiscItem( name="2PL", par=par, est=c(TRUE,TRUE),
                P=P_2PL, prior=c( a="dlnorm"), prior_par1=c(a=0),
                prior_par2=c(a=5) )
item_1PL <- sirt::xxirt_createDiscItem( name="1PL", par=par[2], est=c(TRUE),
                P=P_1PL )
# list of item classes in customItems
customItems <- list( item_1PL,  item_2PL )

#-- create parameter table
itemtype <- rep( "1PL", 12 )
partable <- sirt::xxirt_createParTable(dat, itemtype=itemtype, customItems=customItems)
# privide starting values
partable1 <- sirt::xxirt_modifyParTable( partable, parname="b",
                   value=- stats::qlogis( colMeans(dat) ) )
# equality constraint of parameters and definition of lower bounds
partable1 <- sirt::xxirt_modifyParTable( partable1, item=c("A1","A2"),
                parname="b", parindex=110, lower=-1, value=0)
print(partable1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
