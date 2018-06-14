%% File Name: data.timss.Rd
%% File Version: 0.07

\name{data.timss}
\alias{data.timss}
\docType{data}
\title{
Dataset TIMSS Mathematics
}
\description{
This datasets contains TIMSS mathematics data from 345 students
on 25 items.
}
\usage{data(data.timss)}
\format{
This dataset is a list. \code{data} is the dataset containing
student ID (\code{idstud}), a dummy variable for female (\code{girl})
and student age (\code{age}). The following variables (starting with
\code{M} in the variable name are items.

The format is:

\code{List of 2} \cr
\code{ $ data:'data.frame':} \cr
\code{  ..$ idstud  : num [1:345] 4e+09 4e+09 4e+09 4e+09 4e+09 ...} \cr
\code{  ..$ girl    : int [1:345] 0 0 0 0 0 0 0 0 1 0 ...} \cr
\code{  ..$ age     : num [1:345] 10.5 10 10.25 10.25 9.92 ...} \cr
\code{  ..$ M031286 : int [1:345] 0 0 0 1 1 0 1 0 1 0 ...} \cr
\code{  ..$ M031106 : int [1:345] 0 0 0 1 1 0 1 1 0 0 ...} \cr
\code{  ..$ M031282 : int [1:345] 0 0 0 1 1 0 1 1 0 0 ...} \cr
\code{  ..$ M031227 : int [1:345] 0 0 0 0 1 0 0 0 0 0 ...} \cr
\code{    [...]} \cr
\code{  ..$ M041203 : int [1:345] 0 0 0 1 1 0 0 0 0 1 ...} \cr
\code{ $ item:'data.frame':} \cr
\code{  ..$ item            : Factor w/ 25 levels "M031045","M031068",..: ...} \cr
\code{  ..$ Block           : Factor w/ 2 levels "M01","M02": 1 1 1 1 1 1 ..} \cr
\code{  ..$ Format          : Factor w/ 2 levels "CR","MC": 1 1 1 1 2  ...} \cr
\code{  ..$ Content.Domain  : Factor w/ 3 levels "Data Display",..: 3 3 3 3  ...} \cr
\code{  ..$ Cognitive.Domain: Factor w/ 3 levels "Applying","Knowing",..: 2 3 3 ..} \cr
}




\keyword{datasets}
