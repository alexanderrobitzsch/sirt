%% File Name: data.ml.Rd
%% File Version: 1.07

\name{data.ml}
\alias{data.ml}
\alias{data.ml1}
\alias{data.ml2}

\docType{data}
\title{
Multilevel Datasets
}
\description{
Datasets for conducting multilevel IRT analysis. This dataset is used
in the examples of the function \code{\link{mcmc.2pno.ml}}. 
}
\usage{
data(data.ml1)
data(data.ml2)
}
\format{
\itemize{	
\item  \code{data.ml1}

  A data frame with 2000 student observations in 100 classes on 17 variables.
  The first variable \code{group} contains the class identifier.
  The remaining 16 variables are dichotomous test items. 
  
\code{'data.frame':   2000 obs. of  17 variables:} \cr
\code{ $ group: num  1001 1001 1001 1001 1001 ...} \cr
\code{ $ X1   : num  1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{ $ X2   : num  1 1 1 0 1 1 1 1 1 1 ...} \cr
\code{ $ X3   : num  0 1 1 0 1 0 1 0 1 0 ...} \cr
\code{ $ X4   : num  1 1 1 0 0 1 1 1 1 1 ...} \cr
\code{ $ X5   : num  0 0 0 1 1 1 0 0 1 1 ...} \cr
\code{[...]} \cr
\code{ $ X16  : num  0 0 1 0 0 0 1 0 0 0 ...} \cr

  
\item  \code{data.ml2}

  A data frame with 2000 student observations in 100 classes on 6 variables.
  The first variable \code{group} contains the class identifier.
  The remaining 5 variables are polytomous test items.
  
  \code{'data.frame':   2000 obs. of  6 variables:} \cr
\code{ $ group: num  1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{ $ X1   : num  2 3 4 3 3 3 1 4 4 3 ...} \cr
\code{ $ X2   : num  2 2 4 3 3 2 2 3 4 3 ...} \cr
\code{ $ X3   : num  3 4 5 4 2 3 3 4 4 2 ...} \cr
\code{ $ X4   : num  2 3 3 2 1 3 1 4 4 3 ...} \cr
\code{ $ X5   : num  2 3 3 2 3 3 1 3 2 2 ...} \cr

 }
}




\keyword{datasets}
