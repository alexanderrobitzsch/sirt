%% File Name: automatic.recode.Rd
%% File Version: 0.12

\name{automatic.recode}
\alias{automatic.recode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Automatic Method of Finding Keys in a Dataset with Raw Item Responses
}
\description{
This function calculates keys of a dataset with raw item responses.
It starts with setting the most frequent category of an item to 1.
Then, in each iteration keys are changed such that the highest 
item discrimination is found.
}
\usage{
automatic.recode(data, exclude = NULL, pstart.min = 0.6, allocate = 200, 
    maxiter = 20, progress = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Dataset with raw item responses
}
  \item{exclude}{
Vector with categories to be excluded for searching the key 
}
  \item{pstart.min}{
Minimum probability for an initial solution of keys.
}
  \item{allocate}{
Maximum number of categories per item. This argument is used in
the function \code{tam.ctt3} of the \pkg{TAM} package.
}
  \item{maxiter}{
Maximum number of iterations
}
  \item{progress}{
A logical which indicates if iteration progress should be displayed
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with following entries
\item{item.stat}{Data frame with item name, p value, item discrimination
and the calculated key}
\item{data.scored}{Scored data frame using calculated keys 
in \code{item.stat}
}
\item{categ.stats}{Data frame with statistics for all categories
of all items}
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
#############################################################################
# EXAMPLE 1: data.raw1
#############################################################################
data(data.raw1)

# recode data.raw1 and exclude keys 8 and 9 (missing codes) and
# start with initially setting all categories larger than 50 % frequency to 1
res1 <- sirt::automatic.recode( data.raw1 , exclude=c(8,9) , pstart.min=.50 )
# inspect calculated keys
res1$item.stat

#############################################################################
# EXAMPLE 2: data.timssAusTwn from TAM package
#############################################################################

miceadds::library_install("TAM")
data(data.timssAusTwn,package="TAM")
raw.resp <- data.timssAusTwn[,1:11]
res2 <- sirt::automatic.recode( data=raw.resp )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{R utilities}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
