%% File Name: data.pirlsmissing.Rd
%% File Version: 0.07
%% File Last Change: 2017-01-18 18:08:37

\name{data.pirlsmissing}
\alias{data.pirlsmissing}

\docType{data}

\title{
Dataset from PIRLS Study with Missing Responses
}

\description{
This is a dataset of the PIRLS 2011 study for 4th graders for the reading 
booklet 13 (the 'PIRLS reader') and 4 countries (Austria, Germany, France, Netherlands). 
Missing responses (missing by intention and not reached) are coded by 9.
}

\usage{
data(data.pirlsmissing)
}

\format{
A data frame with 3480 observations on the following 38 variables.

The format is:

\code{'data.frame':   3480 obs. of  38 variables:} \cr
\code{ $ idstud  : int  1000001 1000002 1000003 1000004 1000005 ...} \cr
\code{ $ country : Factor w/ 4 levels "AUT","DEU","FRA",..: 1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{ $ studwgt : num  1.06 1.06 1.06 1.06 1.06 ...} \cr
\code{ $ R31G01M : int  1 1 1 1 1 1 0 1 1 0 ...} \cr
\code{ $ R31G02C : int  0 9 0 1 0 0 0 0 1 0 ...} \cr
\code{ $ R31G03M : int  1 1 1 1 0 1 0 0 1 1 ...} \cr
\code{ [...]} \cr
\code{ $ R31P15C : int  1 9 0 1 0 0 0 0 1 0 ...} \cr
\code{ $ R31P16C : int  0 0 0 0 0 0 0 9 0 1 ...} \cr
      
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}



\examples{
data(data.pirlsmissing)
# inspect missing rates
round( colMeans( data.pirlsmissing==9 ) , 3 )
  ##    idstud  country  studwgt  R31G01M  R31G02C  R31G03M  R31G04C  R31G05M 
  ##     0.000    0.000    0.000    0.009    0.076    0.012    0.203    0.018 
  ##   R31G06M  R31G07M R31G08CZ R31G08CA R31G08CB  R31G09M  R31G10C  R31G11M 
  ##     0.010    0.020    0.189    0.225    0.252    0.019    0.126    0.023 
  ##   R31G12C R31G13CZ R31G13CA R31G13CB R31G13CC  R31G14M  R31P01M  R31P02C 
  ##     0.202    0.170    0.198    0.220    0.223    0.074    0.013    0.039 
  ##   R31P03C  R31P04M  R31P05C  R31P06C  R31P07C  R31P08M  R31P09C  R31P10M 
  ##     0.056    0.012    0.075    0.043    0.074    0.024    0.062    0.025 
  ##   R31P11M  R31P12M  R31P13M  R31P14C  R31P15C  R31P16C 
  ##     0.027    0.030    0.030    0.126    0.130    0.127 
}


\keyword{datasets}
