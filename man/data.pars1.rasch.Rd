%% File Name: data.pars1.rasch.Rd
%% File Version: 0.10
%% File Last Change: 2017-01-18 18:08:37

\name{data.pars1.rasch}
\alias{data.pars1.rasch}
\alias{data.pars1.2pl}

\docType{data}
\title{
Item Parameters for Three Studies Obtained by 1PL and 2PL Estimation
}
\description{
The datasets contain item parameters to be prepared for linking
using the function \code{\link{linking.haberman}}.
}
\usage{
data(data.pars1.rasch)
data(data.pars1.2pl)
}
\format{
\itemize{
\item The format of \code{data.pars1.rasch} is:

\code{'data.frame':   22 obs. of  4 variables:} \cr
\code{ $ study: chr  "study1" "study1" "study1" "study1" ...} \cr
\code{ $ item : Factor w/ 12 levels "M133","M176",..: 1 2 3 4 5 1 6 7 3 8 ...} \cr
\code{ $ a    : num  1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{ $ b    : num  -1.5862 0.40762 1.78031 2.00382 0.00862 ...}

Item slopes \code{a} are fixed to 1 in 1PL estimation. Item difficulties
are denoted by \code{b}. \cr

\item The format of \code{data.pars1.2pl} is:

\code{'data.frame':   22 obs. of  4 variables:} \cr
\code{ $ study: chr  "study1" "study1" "study1" "study1" ...} \cr
\code{ $ item : Factor w/ 12 levels "M133","M176",..: 1 2 3 4 5 1 6 7 3 8 ...} \cr
\code{ $ a    : num  1.238 0.957 1.83 1.927 2.298 ...} \cr
\code{ $ b    : num  -1.16607 0.35844 1.06571 1.17159 0.00792 ...} \cr
}
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}


%data.pars1.rasch
%data.pars1.2pl

%\examples{
%data(data.pars1.rasch)
%## maybe str(data.pars1.rasch) ; plot(data.pars1.rasch) ...
%}
\keyword{datasets}
