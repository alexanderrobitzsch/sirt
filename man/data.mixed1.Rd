%% File Name: data.mixed1.Rd
%% File Version: 0.08

\name{data.mixed1}
\alias{data.mixed1}
\docType{data}

\title{
Dataset with Mixed Dichotomous and Polytomous Item Responses
}

\description{
Dataset with mixed dichotomous and polytomous item responses.
}

\usage{data(data.mixed1)}
\format{
  A data frame with 1000 observations on the following 37 variables.
  
\code{'data.frame':   1000 obs. of  37 variables:} \cr
\code{ $ I01: num  1 1 1 1 1 1 1 0 1 1 ...} \cr
\code{ $ I02: num  1 1 1 1 1 1 1 1 0 1 ...} \cr
\code{ [...]} \cr
\code{ $ I36: num  1 1 1 1 0 0 0 0 1 1 ...} \cr
\code{ $ I37: num  0 1 1 1 0 1 0 0 1 1 ...} \cr
}

%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}

%%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%%}
%%\references{
%%  ~~ possibly secondary sources and usages ~~
%%}



\examples{
data(data.mixed1)
apply( data.mixed1 , 2 , max )
  ##   I01 I02 I03 I04 I05 I06 I07 I08 I09 I10 I11 I12 I13 I14 I15 I16 
  ##     1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
  ##   I17 I18 I19 I20 I21 I22 I23 I24 I25 I26 I27 I28 I29 I30 I31 I32 
  ##     1   1   1   1   4   4   1   1   1   1   1   1   1   1   1   1 
  ##   I33 I34 I35 I36 I37 
  ##     1   1   1   1   1 
}
\keyword{datasets}
