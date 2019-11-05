%% File Name: data.ess2005.Rd
%% File Version: 0.11

\name{data.ess2005}
\alias{data.ess2005}
\docType{data}
\title{
Dataset European Social Survey 2005
}
\description{
This dataset contains item loadings \eqn{\lambda} and intercepts \eqn{\nu}
for 26 countries for the European Social Survey (ESS 2005;
see Asparouhov & Muthen, 2014).
}

\usage{
data(data.ess2005)
}

\format{
The format of the dataset is:

\code{List of 2} \cr
\code{ $ lambda: num [1:26, 1:4] 0.688 0.721 0.72 0.687 0.625 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 2} \cr
\code{  .. ..$ : NULL} \cr
\code{  .. ..$ : chr [1:4] "ipfrule" "ipmodst" "ipbhprp" "imptrad"} \cr
\code{ $ nu    : num [1:26, 1:4] 3.26 2.52 3.41 2.84 2.79 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 2} \cr
\code{  .. ..$ : NULL} \cr
\code{  .. ..$ : chr [1:4] "ipfrule" "ipmodst" "ipbhprp" "imptrad"} \cr

}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}

% \section{Example Index}{
% \code{\link{invariance.alignment}} (Example 3)
% }

\references{
Asparouhov, T., & Muthen, B. (2014). Multiple-group factor analysis alignment.
\emph{Structural Equation Modeling, 21}(4), 1-14.
\doi{10.1080/10705511.2014.919210}
}
%\examples{
%data(data.ess2005)
%## maybe str(data.ess2005) ; plot(data.ess2005) ...
%}


%% \keyword{datasets}
