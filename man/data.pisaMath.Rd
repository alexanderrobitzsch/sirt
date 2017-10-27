%% File Name: data.pisaMath.Rd
%% File Version: 0.08

\name{data.pisaMath}
\alias{data.pisaMath}
\docType{data}
\title{
Dataset PISA Mathematics
}
\description{
This is an example PISA dataset of mathematics items. The
dataset contains 565 students on 11 items.
}
\usage{
data(data.pisaMath)
}


\format{
The dataset is a list. The list element \code{data}
contains the dataset with the demographical variables 
student ID (\code{idstud}), school ID (\code{idschool}), a dummy variable for female 
students (\code{female}), socioeconomic status (\code{hisei})
and migrational background (\code{migra}). The remaining 
variables (starting with \code{M} in the name) are
the mathematics items. \cr
The item metadata are included in the list element
\code{item} which contains item name (\code{item}) and the
testlet label (\code{testlet}). An item not included
in a testlet is indicated by \code{NA}.

The format is:

\code{List of 2} \cr
\code{ $ data:'data.frame':} \cr
\code{  ..$ idstud  : num [1:565] 9e+10 9e+10 9e+10 9e+10 9e+10 ...} \cr
\code{  ..$ idschool: int [1:565] 900015 900015 900015 900015  ...} \cr
\code{  ..$ female  : int [1:565] 0 0 0 0 0 0 0 0 0 0 ...} \cr
\code{  ..$ hisei   : num [1:565] -1.16 -1.099 -1.588 -0.365 -1.588 ...} \cr
\code{  ..$ migra   : int [1:565] 0 0 0 0 0 0 0 0 0 1 ...} \cr
\code{  ..$ M192Q01 : int [1:565] 1 0 1 1 1 1 1 0 0 0 ...} \cr
\code{  ..$ M406Q01 : int [1:565] 1 1 1 0 1 0 0 0 1 0 ...} \cr
\code{  ..$ M406Q02 : int [1:565] 1 0 0 0 1 0 0 0 1 0 ...} \cr
\code{  ..$ M423Q01 : int [1:565] 0 1 0 1 1 1 1 1 1 0 ...} \cr
\code{  ..$ M496Q01 : int [1:565] 1 0 0 0 0 0 0 0 1 0 ...} \cr
\code{  ..$ M496Q02 : int [1:565] 1 0 0 1 0 1 0 1 1 0 ...} \cr
\code{  ..$ M564Q01 : int [1:565] 1 1 1 1 1 1 0 0 1 0 ...} \cr
\code{  ..$ M564Q02 : int [1:565] 1 0 1 1 1 0 0 0 0 0 ...} \cr
\code{  ..$ M571Q01 : int [1:565] 1 0 0 0 1 0 0 0 0 0 ...} \cr
\code{  ..$ M603Q01 : int [1:565] 1 0 0 0 1 0 0 0 0 0 ...} \cr
\code{  ..$ M603Q02 : int [1:565] 1 0 0 0 1 0 0 0 1 0 ...} \cr
\code{ $ item:'data.frame':} \cr
\code{  ..$ item   : Factor w/ 11 levels "M192Q01","M406Q01",..: 1 2 3 4  ...} \cr
\code{  ..$ testlet: chr [1:11] NA "M406" "M406" NA ...} \cr
}



\keyword{datasets}
