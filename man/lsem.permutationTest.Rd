%% File Name: lsem.permutationTest.Rd
%% File Version: 0.18

\name{lsem.permutationTest}
\alias{lsem.permutationTest}
\alias{summary.lsem.permutationTest}
\alias{plot.lsem.permutationTest}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Permutation Test for a Local Structural Equation Model
}

\description{
Performs a permutation test for testing the hypothesis that model
parameter are independent of a moderator variable (see Hildebrandt,
Wilhelm, & Robitzsch, 2009).
}


\usage{
lsem.permutationTest(lsem.object, B = 1000, residualize = TRUE, verbose = TRUE)

\method{summary}{lsem.permutationTest}(object, file=NULL, digits=3, ...)

\method{plot}{lsem.permutationTest}(x, type = "global", stattype = "SD", 
    parindex = NULL, sig_add = TRUE, sig_level = 0.05, sig_pch=17, nonsig_pch=2, 
    sig_cex = 1, sig_lab = "p value",  stat_lab = "Test statistic", 
    moderator_lab = NULL, digits = 3, title = NULL, parlabels = NULL, 
    ask = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lsem.object}{
Fitted object of class \code{lsem} with \code{\link{lsem.estimate}}
}
  \item{B}{
Number of permutation samples
}
  \item{residualize}{
Optional logical indicating whether residualization of the moderator
should be performed for each permutation sample. 
}
  \item{verbose}{
Optional logical printing information about computation progress.
}
\item{object}{Object of class \code{lsem}}
\item{file}{A file name in which the summary output will be written.}
\item{digits}{Number of digits.}
\item{\dots}{Further arguments to be passed.}
  \item{x}{
Object of class \code{lsem}
}
  \item{type}{
Type of the statistic to be plotted. If \code{type="global"}, a global
test will be displayed. If \code{type="pointwise"} for each value at the 
focal point (defined in \code{moderator.grid}) are calculated.
}
  \item{stattype}{
Type of test statistics. Can be \code{MAD} (mean absolute deviation), 
\code{SD} (standard deviation) or \code{lin_slo} (linear slope). 
}
  \item{parindex}{
Vector of indices of selected parameters. 
}
  \item{sig_add}{
Logical indicating whether significance values (p values) should
be displayed.
}
  \item{sig_level}{
Significance level.
}
  \item{sig_pch}{
Point symbol for significant values.
}
  \item{nonsig_pch}{
Point symbol for non-significant values.
}
\item{sig_cex}{Point size for graphic displaying p values}
  \item{sig_lab}{
Label for significance value (p value).
}
\item{stat_lab}{Label of y axis for graphic with pointwise test statistic}
  \item{moderator_lab}{
Label of the moderator.
}
  \item{title}{
Title of the plot. Can be a vector.
}
  \item{parlabels}{
Labels of the parameters. Can be a vector.
}
\item{ask}{A logical which asks for changing the graphic for each parameter.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
List with following entries
\item{teststat}{Data frame with global test statistics. The statistics
are \code{SD}, \code{MAD} and \code{lin_slo} with their corresponding
p values.}
\item{parameters_pointwise_test}{Data frame with pointwise test statistics.}
\item{parameters}{Original parameters.}
\item{parameters}{Parameters in permutation samples.}
\item{parameters_summary}{Original parameter summary.}
\item{parameters_summary_M}{Mean of each parameter in permutation sample.}
\item{parameters_summary_SD}{Standard deviation (SD) statistic in permutation slope.}
\item{parameters_summary_MAD}{Mean absolute deviation (MAD)
	statistic in permutation sample.}
\item{parameters_summary_MAD}{Linear slope parameter in permutation sample.}
}

%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...


\references{
Hildebrandt, A., Wilhelm, O., & Robitzsch, A. (2009). Complementary and 
competing factor analytic approaches for the investigation of measurement invariance. 
\emph{Review of Psychology}, \bold{16}, 87-102.
}

\author{
Alexander Robitzsch, Oliver Luedtke, Andrea Hildebrandt
}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For Examples see \code{\link{lsem.estimate}}.
}

%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{LSEM}
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
