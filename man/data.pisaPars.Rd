\name{data.pisaPars}
\alias{data.pisaPars}
\docType{data}
\title{
Item Parameters from Two PISA Studies
}
\description{
This data frame contains item parameters from two PISA studies. 
Because the Rasch model is used, only item difficulties are considered.
}
\usage{
data(data.pisaPars)
}

\format{
  A data frame with 25 observations on the following 4 variables. 
  \describe{
    \item{\code{item}}{Item names}
    \item{\code{testlet}}{Items are arranged in corresonding testlets. These
		names are located in this column.}
    \item{\code{study1}}{Item difficulties of study 1}
    \item{\code{study2}}{Item difficulties of study 2}
		}
}


\keyword{datasets}
