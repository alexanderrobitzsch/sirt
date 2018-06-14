%% File Name: categorize.Rd
%% File Version: 0.11

\name{categorize}
\alias{categorize}
\alias{decategorize}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Categorize and Decategorize Variables in a Data Frame
}

\description{
The function \code{categorize} defines categories for variables in
a data frame, starting with a user-defined index (e.g. 0 or 1).
Continuous variables can be categorized by defining categories by
discretizing the variables in different quantile groups.

The function \code{decategorize} does the reverse operation.
}


\usage{
categorize(dat, categorical = NULL, quant=NULL , lowest = 0)

decategorize(dat, categ_design = NULL)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{dat}{
Data frame
}
  \item{categorical}{
Vector with variable names which should be converted into categories,
beginning with integer \code{lowest}
}
\item{quant}{Vector with number of classes for each variables.
    Variables are categorized among quantiles. The vector must
have names containing variable names.}
  \item{lowest}{
Lowest category index. Default is 0.
}
\item{categ_design}{Data frame containing informations about
categorization which is the output of \code{categorize}.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
For \code{categorize}, it is a list with entries
  \item{data}{Converted data frame}
  \item{categ_design}{Data frame containing some informations
    about categorization}

For \code{decategorize} it is a data frame.
}

%\references{
%% ~put references to the literature/web site here ~
%}
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
library(mice)
library(miceadds)

#############################################################################
# EXAMPLE 1: Categorize questionnaire data
#############################################################################

data(data.smallscale , package="miceadds")
dat <- data.smallscale

# (0) select dataset
dat <- dat[ , 9:20 ]
summary(dat)
categorical <- colnames(dat)[2:6]

# (1) categorize data
res <- sirt::categorize( dat , categorical=categorical )

# (2) multiple imputation using the mice package
dat2 <- res$data
VV <- ncol(dat2)
impMethod <- rep( "sample" , VV )    # define random sampling imputation method
names(impMethod) <- colnames(dat2)
imp <- mice::mice( as.matrix(dat2) , impMethod = impMethod , maxit=1 , m=1 )
dat3 <- mice::complete(imp,action=1)

# (3) decategorize dataset
dat3a <- sirt::decategorize( dat3 , categ_design = res$categ_design )

#############################################################################
# EXAMPLE 2: Categorize ordinal and continuous data
#############################################################################

data(data.ma01,package="miceadds")
dat <- data.ma01
summary(dat[,-c(1:2)] )

# define variables to be categorized
categorical <- c("books" , "paredu" )
# define quantiles
quant <-  c(6,5,11)
names(quant) <- c("math" , "read" , "hisei")

# categorize data
res <- sirt::categorize( dat , categorical = categorical , quant=quant)
str(res)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Utilities}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
