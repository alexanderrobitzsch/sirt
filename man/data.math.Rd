%% File Name: data.math.Rd
%% File Version: 0.12

\name{data.math}
\alias{data.math}
\docType{data}
\title{
Dataset Mathematics
}
\description{
This is an example Mathematics dataset. The
dataset contains 664 students on 30 items.
}
\usage{data(data.math)}
\format{
The dataset is a list. The list element \code{data}
contains the dataset with the demographical variables
student ID (\code{idstud}) and a dummy variable
for female students (\code{female}). The remaining
variables (starting with \code{M} in the name) are
the mathematics items. \cr
The item metadata are included in the list element
\code{item} which contains item name (\code{item}) and the
testlet label (\code{testlet}). An item not included
in a testlet is indicated by \code{NA}.
Each item is allocated to one and only competence domain (\code{domain}).
\cr

The format is:

\code{List of 2} \cr
\code{ $ data:'data.frame':} \cr
\code{  ..$ idstud: int [1:664] 1001 1002 1003 ...} \cr
\code{  ..$ female: int [1:664] 1 1 0 0 1 1 1 0 0 1 ...} \cr
\code{  ..$ MA1   : int [1:664] 1 1 1 0 0 1 1 1 1 1 ...} \cr
\code{  ..$ MA2   : int [1:664] 1 1 1 1 1 0 0 0 0 1 ...} \cr
\code{  ..$ MA3   : int [1:664] 1 1 0 0 0 0 0 1 0 0 ...} \cr
\code{  ..$ MA4   : int [1:664] 0 1 1 1 0 0 1 0 0 0 ...} \cr
\code{  ..$ MB1   : int [1:664] 0 1 0 1 0 0 0 0 0 1 ...} \cr
\code{  ..$ MB2   : int [1:664] 1 1 1 1 0 1 0 1 0 0 ...} \cr
\code{  ..$ MB3   : int [1:664] 1 1 1 1 0 0 0 1 0 1 ...} \cr
\code{  [...]} \cr
\code{  ..$ MH3   : int [1:664] 1 1 0 1 0 0 1 0 1 0 ...} \cr
\code{  ..$ MH4   : int [1:664] 0 1 1 1 0 0 0 0 1 0 ...} \cr
\code{  ..$ MI1   : int [1:664] 1 1 0 1 0 1 0 0 1 0 ...} \cr
\code{  ..$ MI2   : int [1:664] 1 1 0 0 0 1 1 0 1 1 ...} \cr
\code{  ..$ MI3   : int [1:664] 0 1 0 1 0 0 0 0 0 0 ...} \cr
\code{ $ item:'data.frame':} \cr
\code{  ..$ item     : Factor w/ 30 levels "MA1","MA2","MA3",..: 1 2 3 4 5 ...} \cr
\code{  ..$ testlet  : Factor w/ 9 levels "","MA","MB","MC",..: 2 2 2 2 3 3  ...} \cr
\code{  ..$ domain   : Factor w/ 3 levels "arithmetic","geometry",..: 1 1 1  ...} \cr
\code{  ..$ subdomain: Factor w/ 9 levels "","addition",..: 2 2 2 2 7 7  ...} \cr
}





%% \keyword{datasets}
