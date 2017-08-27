%% File Name: data.ratings1.Rd
%% File Version: 0.17
%% File Last Change: 2017-01-18 18:08:37

\name{data.ratings}
\alias{data.ratings}
\alias{data.ratings1}
\alias{data.ratings2}
\alias{data.ratings3}


\docType{data}
\title{
Rating Datasets
}
\description{
Some rating datasets.
}
\usage{
data(data.ratings1)
data(data.ratings2)
data(data.ratings3)
}



\format{
	
\itemize{	
\item Dataset \code{data.ratings1}: \cr

Data frame with 274 observations containing 5 criteria (\code{k1}, ..., \code{k5}),
135 students and 7 raters.
  
\code{'data.frame':   274 obs. of  7 variables:} \cr
\code{ $ idstud: int  100020106 100020106 100070101 100070101 100100109  ...} \cr
\code{ $ rater : Factor w/ 16 levels "db01","db02",..: 3 15 5 10 2 1 5 4 1 5 ...} \cr
\code{ $ k1    : int  1 1 0 1 2 0 1 3 0 0 ...} \cr
\code{ $ k2    : int  1 1 1 1 1 0 0 3 0 0 ...} \cr
\code{ $ k3    : int  1 1 1 1 2 0 0 3 1 0 ...} \cr
\code{ $ k4    : int  1 1 1 2 1 0 0 2 0 1 ...} \cr
\code{ $ k5    : int  2 2 1 2 0 1 0 3 1 0 ...} \cr

Data from a 2009 Austrian survey of national educational
standards for 8th graders in German language writing.
Variables \code{k1} to \code{k5} denote several rating
criteria of writing compentency. \cr

\item Dataset \code{data.ratings2}: \cr

Data frame with 615 observations containing 5 criteria (\code{k1}, ..., \code{k5}),
178 students and 16 raters.

\code{'data.frame':   615 obs. of  7 variables:} \cr
\code{ $ idstud: num  1001 1001 1002 1002 1003 ...} \cr
\code{ $ rater : chr  "R03" "R15" "R05" "R10" ...} \cr
\code{ $ k1    : int  1 1 0 1 2 0 1 3 3 0 ...} \cr
\code{ $ k2    : int  1 1 1 1 1 0 0 3 3 0 ...} \cr
\code{ $ k3    : int  1 1 1 1 2 0 0 3 3 1 ...} \cr
\code{ $ k4    : int  1 1 1 2 1 0 0 2 2 0 ...} \cr
\code{ $ k5    : int  2 2 1 2 0 1 0 3 2 1 ...} \cr

\item Dataset \code{data.ratings3}: \cr

Data frame with 3169 observations containing 4 criteria (\code{crit2}, ..., \code{crit6}),
561 students and 52 raters.

\code{'data.frame':   3169 obs. of  6 variables:} \cr
\code{ $ idstud: num  10001 10001 10002 10002 10003 ...} \cr
\code{ $ rater : num  840 838 842 808 830 845 813 849 809 802 ...} \cr
\code{ $ crit2 : int  1 3 3 1 2 2 2 2 3 3 ...} \cr
\code{ $ crit3 : int  2 2 2 2 2 2 2 2 3 3 ...} \cr
\code{ $ crit4 : int  1 2 2 2 1 1 1 2 2 2 ...} \cr
\code{ $ crit6 : num  4 4 4 3 4 4 4 4 4 4 ...} \cr


}
}

%\source{
%xxx
%}






%\examples{
%# include an example of rating data here	
%}	

\keyword{datasets}
