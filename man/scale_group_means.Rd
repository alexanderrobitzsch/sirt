%% File Name: scale_group_means.Rd
%% File Version: 0.03

\name{scale_group_means}
\alias{scale_group_means}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Scaling of Group Means and Standard Deviations
}

\description{
Scales a vector of means and standard deviations containing group values.
}


\usage{
scale_group_means(M, SD, probs=NULL)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{M}{Vector of means}
\item{SD}{Vector of standard deviations}
\item{probs}{Optional vector containing probabilities}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
For \code{categorize}, it is a list with entries
  \item{data}{Converted data frame}
  \item{categ_design}{Data frame containing some informations
    about categorization}

For \code{decategorize} it is a data frame.
}

%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%Alexander Robitzsch
%}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~


%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
#############################################################################
# EXAMPLE 1: Toy example
#############################################################################

M <- c(-.03, .18, -.23, -.15, .29)
SD <- c(.97, 1.13, .77, 1.05, 1.17)
sirt::scale_group_means(M=M, SD=SD)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Utilities}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
