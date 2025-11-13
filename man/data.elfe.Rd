%% File Name: data.elfe.Rd
%% File Version: 0.031

\name{data.elfe}
\alias{data.elfe}
\docType{data}
\title{
ELFE Dataset
}
\description{
This is synthetic dataset of the ELFE reading comprehension test used
in Robitzsch et al. (2011).
}
\usage{
data(data.elfe)
}

\format{
A data frame with 1855 observations on the following 59 variables. 
Included are item responses to 19 items administered at three time
points.

The format is

\code{'data.frame':   1855 obs. of  59 variables:} \cr
\code{ $ index : int  1 2 3 4 5 6 7 8 9 10 ...} \cr
\code{ $ idstud: int  100001 100002 100003 100004 100005 100006 100007 100008 100009 100010 ...} \cr
\code{ $ T1I02 : int  1 0 1 1 1 0 1 1 0 0 ...} \cr
\code{[...]} \cr
\code{ $ T1I20 : int  1 0 0 0 0 0 0 0 0 0 ...} \cr
\code{ $ T2I02 : int  1 1 1 1 1 1 1 1 0 1 ...} \cr
\code{[...]} \cr
\code{ $ T2I20 : int  0 0 1 1 0 0 0 0 0 0 ...} \cr
\code{ $ T3I02 : int  1 1 1 1 0 1 1 1 1 1 ...} \cr
\code{[...]} \cr
\code{ $ T3I20 : int  0 0 1 0 0 0 0 1 0 0 ...} \cr
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{
Robitzsch, A., Doerfler, T., Pfost, M., & Artelt, C. (2011). Die Bedeutung der Itemauswahl 
und der Modellwahl fuer die laengsschnittliche Erfassung von Kompetenzen: 
Lesekompetenzentwicklung in der Primarstufe. 
\emph{Zeitschrift fuer Entwicklungspsychologie und paedagogische Psychologie, 43}, 
213-227. \doi{10.1026/0049-8637/a000052}
}

%\examples{
%data(data.elfe)
%## maybe str(data.elfe) ; plot(data.elfe) ...
%}
%\keyword{datasets}
