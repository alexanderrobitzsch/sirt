\name{data.pisaRead}
\alias{data.pisaRead}
\docType{data}
\title{
Dataset PISA Reading
}
\description{
This is an example PISA dataset of reading items. The
dataset contains 623 students on 12 items.
}
\usage{data(data.pisaRead)}
\format{
The dataset is a list. The list element \code{data}
contains the dataset with the demographical variables 
student ID (\code{idstud}), school ID (\code{idschool}), a dummy variable for female 
students (\code{female}), socioeconomic status (\code{hisei})
and migrational background (\code{migra}). The remaining 
variables (starting with \code{R} in the name) are
the reading items. \cr
The item metadata are included in the list element
\code{item} which contains item name (\code{item}),
testlet label (\code{testlet}), item format (\code{ItemFormat}),
text type (\code{TextType}) and text aspect (\code{Aspect}). 

The format is:

\code{List of 2} \cr
\code{ $ data:'data.frame':} \cr
\code{  ..$ idstud  : num [1:623] 9e+10 9e+10 9e+10 9e+10 9e+10 ...} \cr
\code{  ..$ idschool: int [1:623] 900003 900003 900003 900003 ...} \cr
\code{  ..$ female  : int [1:623] 1 0 1 0 0 0 1 0 1 0 ...} \cr
\code{  ..$ hisei   : num [1:623] -1.16 -0.671 1.286 0.185 1.225 ...} \cr
\code{  ..$ migra   : int [1:623] 0 0 0 0 0 0 0 0 0 0 ...} \cr
\code{  ..$ R432Q01 : int [1:623] 1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{  ..$ R432Q05 : int [1:623] 1 1 1 1 1 0 1 1 1 0 ...} \cr
\code{  ..$ R432Q06 : int [1:623] 0 0 0 0 0 0 0 0 0 0 ...} \cr
\code{  ..$ R456Q01 : int [1:623] 1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{  ..$ R456Q02 : int [1:623] 1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{  ..$ R456Q06 : int [1:623] 1 1 1 1 1 1 0 0 1 1 ...} \cr
\code{  ..$ R460Q01 : int [1:623] 1 1 0 0 0 0 0 1 1 1 ...} \cr
\code{  ..$ R460Q05 : int [1:623] 1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{  ..$ R460Q06 : int [1:623] 0 1 1 1 1 1 0 0 1 1 ...} \cr
\code{  ..$ R466Q02 : int [1:623] 0 1 0 1 1 0 1 0 0 1 ...} \cr
\code{  ..$ R466Q03 : int [1:623] 0 0 0 1 0 0 0 1 0 1 ...} \cr
\code{  ..$ R466Q06 : int [1:623] 0 1 1 1 1 1 0 1 1 1 ...} \cr
\code{ $ item:'data.frame':} \cr
\code{  ..$ item      : Factor w/ 12 levels "R432Q01","R432Q05",..: 1 2 3 4  ...} \cr
\code{  ..$ testlet   : Factor w/ 4 levels "R432","R456",..: 1 1 1 2  ...} \cr
\code{  ..$ ItemFormat: Factor w/ 2 levels "CR","MC": 1 1 2 2 1 1 1 2 2 1 ...} \cr
\code{  ..$ TextType  : Factor w/ 3 levels "Argumentation",..: 1 1 1 3  ...} \cr
\code{  ..$ Aspect    : Factor w/ 3 levels "Access_and_retrieve",..: 2 3 2 1 ...} \cr
}




\keyword{datasets}
