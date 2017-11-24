%% File Name: data.befki.Rd
%% File Version: 0.07

\name{data.befki}
\alias{data.befki}
\alias{data.befki_resp}

\docType{data}

\title{
BEFKI Dataset (Schroeders, Schipolowski, & Wilhelm, 2015)
}

\description{
The synthetic dataset is based on the standardization sample of the Berlin
Test of Fluid and Crystallized Intelligence (BEFKI, Wilhelm, Schroeders, &
Schipolowski, 2014). The underlying sample consists of N = 11,756 students
from all German federal states (except for the smallest one) and all school
types of the general educational system attending Grades 5 to 12. A
detailed description of the study, the sample, and the measure is given in
Schroeders, Schipolowski, and Wilhelm (2015).
}


\usage{
data(data.befki)
data(data.befki_resp)
}


\format{
\itemize{	
\item The dataset \code{data.befki} contains 11756 students, nested within
581 classes.

\code{'data.frame':   11756 obs. of  12 variables:} \cr
\code{ $ idclass: int  1276 1276 1276 1276 1276 1276 1276 1276 1276 1276 ...} \cr
\code{ $ idstud : int  127601 127602 127603 127604 127605 127606 127607 127608 127609 127610 ...} \cr
\code{ $ grade  : int  5 5 5 5 5 5 5 5 5 5 ...} \cr
\code{ $ gym    : int  0 0 0 0 0 0 0 0 0 0 ...} \cr
\code{ $ female : int  0 1 0 0 0 0 1 0 0 0 ...} \cr
\code{ $ age    : num  12.2 11.8 11.5 10.8 10.9 ...} \cr
\code{ $ sci    : num  -3.14 -3.44 -2.62 -2.16 -1.01 -1.91 -1.01 -4.13 -2.16 -3.44 ...} \cr
\code{ $ hum    : num  -1.71 -1.29 -2.29 -2.48 -0.65 -0.92 -1.71 -2.31 -1.99 -2.48 ...} \cr
\code{ $ soc    : num  -2.87 -3.35 -3.81 -2.35 -1.32 -1.11 -1.68 -2.96 -2.69 -3.35 ...} \cr
\code{ $ gfv    : num  -2.25 -2.19 -2.25 -1.17 -2.19 -3.05 -1.7 -2.19 -3.05 -1.7 ...} \cr
\code{ $ gfn    : num  -2.2 -1.85 -1.85 -1.85 -1.85 -0.27 -1.37 -2.58 -1.85 -3.13 ...} \cr
\code{ $ gff    : num  -0.91 -0.43 -1.17 -1.45 -0.61 -1.78 -1.17 -1.78 -1.78 -3.87 ...} \cr

\item The dataset \code{data.befki_resp} contains response indicators for observed data points
in the dataset \code{data.befki}.

\code{ num [1:11756, 1:12] 1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{ - attr(*, "dimnames")=List of 2} \cr
\code{  ..$ : NULL} \cr
\code{  ..$ : chr [1:12] "idclass" "idstud" "grade" "gym" ...} \cr

}


}

\details{
The procedure for generating this dataset is based on a factorization of the
joint distribution. All variables are simulated from unidimensional conditional
parametric regression models including several interaction and quadratic terms. 
The multilevel structure is approximated by including
cluster means as predictors in the regression models. 
}

\source{
Synthetic dataset
}


\references{
Schroeders, U., Schipolowski, S., & Wilhelm, O. (2015). Age-related
changes in the mean and covariance structure of fluid and crystallized
intelligence in childhood and adolescence. \emph{Intelligence}, \bold{48}, 15-29.
DOI 10.1016/j.intell.2014.10.006

Wilhelm, O., Schroeders, U., & Schipolowski, S. (2014). \emph{Berliner Test
zur Erfassung fluider und kristalliner Intelligenz fuer die 8. bis 10.
Jahrgangsstufe} [Berlin test of fluid and crystallized intelligence for
grades 8-10]. Goettingen: Hogrefe.
}


%\examples{
%data(data.befki)
%}
%\keyword{datasets}
