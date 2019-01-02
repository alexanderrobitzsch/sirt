%% File Name: prior_model_parse.Rd
%% File Version: 0.10

\name{prior_model_parse}
\alias{prior_model_parse}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Parsing a Prior Model
}
\description{
Parses a string specifying a prior model which is needed
for the \code{prior} argument in \code{LAM::amh}
}
\usage{
prior_model_parse(prior_model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prior_model}{
String specifying the prior conforming to \R syntax.
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
List with specified prior distributions for parameters
as needed for the \code{prior} argument in \code{LAM::amh}
}
%

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
\code{LAM::amh}
}

\examples{
#############################################################################
# EXAMPLE 1: Toy example prior distributions
#############################################################################

#*** define prior model as a string
prior_model <- "
  # prior distributions means
  mu1 ~ dnorm( NA, mean=0, sd=1 )
  mu2 ~ dnorm(NA)       # mean T2 and T3
  # prior distribution standard deviation
  sig1 ~ dunif(NA,0, max=10)
      "

#*** convert priors into a list
res <- sirt::prior_model_parse( prior_model )
str(res)
  ##  List of 3
  ##   $ mu1 :List of 2
  ##    ..$ : chr "dnorm"
  ##    ..$ :List of 3
  ##    .. ..$ NA  : num NA
  ##    .. ..$ mean: num 0
  ##    .. ..$ sd  : num 1
  ##   $ mu2 :List of 2
  ##    ..$ : chr "dnorm"
  ##    ..$ :List of 1
  ##    .. ..$ : num NA
  ##   $ sig1:List of 2
  ##    ..$ : chr "dunif"
  ##    ..$ :List of 3
  ##    .. ..$ NA : num NA
  ##    .. ..$ NA : num 0
  ##    .. ..$ max: num 10
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
