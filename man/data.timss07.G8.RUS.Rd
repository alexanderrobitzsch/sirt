%% File Name: data.timss07.G8.RUS.Rd
%% File Version: 0.03

\name{data.timss07.G8.RUS}
\alias{data.timss07.G8.RUS}
\docType{data}

\title{
TIMSS 2007 Grade 8 Mathematics and Science Russia
}

\description{
This TIMSS 2007 dataset contains item responses of 4472 eigth grade
Russian students in Mathematics and Science.
}


\usage{
data(data.timss07.G8.RUS)
}


\format{
The datasets contains raw responses (\code{raw}), scored responses
(\code{scored}) and item informations (\code{iteminfo}).
	
The format of the dataset is:

\code{List of 3} \cr
\code{ $ raw     :'data.frame':} \cr
\code{  ..$ idstud  : num [1:4472] 3010101 3010102 3010104 3010105 3010106 ...} \cr
\code{  ..$ M022043 : atomic [1:4472] NA 1 4 NA NA NA NA NA NA NA ...} \cr
\code{  .. ..- attr(*, "value.labels")= Named num [1:7] 9 6 5 4 3 2 1} \cr
\code{  .. .. ..- attr(*, "names")= chr [1:7] "OMITTED" "NOT REACHED" "E" "D*" ...} \cr
\code{  [...]} \cr
\code{  ..$ M032698 : atomic [1:4472] NA NA NA NA NA NA NA 2 1 NA ...} \cr
\code{  .. ..- attr(*, "value.labels")= Named num [1:6] 9 6 4 3 2 1} \cr
\code{  .. .. ..- attr(*, "names")= chr [1:6] "OMITTED" "NOT REACHED" "D" "C" ...} \cr
\code{  ..$ M032097 : atomic [1:4472] NA NA NA NA NA NA NA 2 3 NA ...} \cr
\code{  .. ..- attr(*, "value.labels")= Named num [1:6] 9 6 4 3 2 1} \cr
\code{  .. .. ..- attr(*, "names")= chr [1:6] "OMITTED" "NOT REACHED" "D" "C*" ...} \cr
\code{  .. [list output truncated]} \cr
\code{ $ scored  : num [1:4472, 1:443] 3010101 3010102 3010104 3010105 3010106 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 2} \cr
\code{  .. ..$ : NULL} \cr
\code{  .. ..$ : chr [1:443] "idstud" "M022043" "M022046" "M022049" ...} \cr
\code{ $ iteminfo:'data.frame':} \cr
\code{  ..$ item      : Factor w/ 442 levels "M022043","M022046",..: 1 2 3 4 5 6 21 7 8 17 ...} \cr
\code{  ..$ content   : Factor w/ 8 levels "Algebra","Biology",..: 7 7 6 1 6 7 4 6 7 7 ...} \cr
\code{  ..$ topic     : Factor w/ 49 levels "Algebraic Expression",..: 32 32 41 29  ...} \cr
\code{  ..$ cognitive : Factor w/ 3 levels "Applying","Knowing",..: 2 1 3 2 1 1 1 1 2 1 ...} \cr
\code{  ..$ item.type : Factor w/ 2 levels "CR","MC": 2 1 2 2 1 2 2 2 2 1 ...} \cr
\code{  ..$ N.options : Factor w/ 4 levels "-"," -","4","5": 4 1 3 4 1 4 4 4 3 1 ...} \cr
\code{  ..$ key       : Factor w/ 7 levels "-"," -","A","B",..: 6 1 6 7 1 5 5 4 6 1 ...} \cr
\code{  ..$ max.points: int [1:442] 1 1 1 1 1 1 1 1 1 2 ...} \cr
\code{  ..$ item.label: Factor w/ 432 levels "1 teacher for every 12 students ",..: 58 351 ...} \cr
  
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
TIMSS 2007 8th Grade, Russian Sample
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}



%\section{Example Index}{
%\code{\link{automatic.recode}} (Example 1),
%}


%\examples{
%data(data.timss07.G8.RUS)
%## maybe str(data.timss07.G8.RUS) ; plot(data.timss07.G8.RUS) ...
%}


\keyword{datasets}
