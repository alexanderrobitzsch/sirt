\name{data.activity.itempars}
\alias{data.activity.itempars}
\docType{data}
\title{
Item Parameters Cultural Activities
}
\description{
List with item parameters for cultural activities
of Austrian students for 9 Austrian countries.
}
\usage{
data(data.activity.itempars)
}

\format{
The format is a list with number of students per group
(\code{N}), item loadings (\code{lambda}) and
item intercepts (\code{nu}):

\code{List of 3} \cr
\code{ $ N     : 'table' int [1:9(1d)] 2580 5279 15131 14692 5525 11005 7080 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 1} \cr
\code{  .. ..$ : chr [1:9] "1" "2" "3" "4" ...} \cr
\code{ $ lambda: num [1:9, 1:5] 0.423 0.485 0.455 0.437 0.502 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 2} \cr
\code{  .. ..$ : chr [1:9] "country1" "country2" "country3" "country4" ...} \cr
\code{  .. ..$ : chr [1:5] "act1" "act2" "act3" "act4" ...} \cr
\code{ $ nu    : num [1:9, 1:5] 1.65 1.53 1.7 1.59 1.7 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 2} \cr
\code{  .. ..$ : chr [1:9] "country1" "country2" "country3" "country4" ...} \cr
\code{  .. ..$ : chr [1:5] "act1" "act2" "act3" "act4" ...} \cr
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%% }
%% \source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%%}
%%\references{
%%  ~~ possibly secondary sources and usages ~~
%% }



%%\examples{
%%data(data.activity.itempars)
%%## maybe str(data.activity.itempars) ; plot(data.activity.itempars) ...
%%}
\keyword{datasets}
