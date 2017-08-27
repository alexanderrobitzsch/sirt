%% File Name: parmsummary_extend.Rd
%% File Version: 0.05
%% File Last Change: 2017-07-11 16:35:55

\name{parmsummary_extend}
\alias{parmsummary_extend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Includes Confindence Interval in Parameter Summary Table 
}
\description{
Includes confindence interval in parameter summary table.
}
\usage{
parmsummary_extend(dfr, level = .95 , est_label = "est", se_label = "se", 
      df_label = "df")
}
%- maybe also 'usage' for other objects documented here.


\arguments{
  \item{dfr}{
Data frame containing parameter summary
}
  \item{level}{
Significance level
}
  \item{est_label}{
Label for parameter estimate
}
  \item{se_label}{
Label for standard error
}
  \item{df_label}{
Label for degrees of freedom
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
Extended parameter summary table
}

%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[stats:confint]{stats::confint}}
}


\examples{
#############################################################################
## EXAMPLE 1: Toy example parameter summary table 
#############################################################################

dfr <- data.frame( "parm" = c("b0" , "b1" ) , "est" = c(0.1 , 1.3 ) , 
                "se" = c(.21 , .32) )
print( sirt::parmsummary_extend(dfr) , digits= 4 )
  ##    parm est   se      t         p lower95 upper95
  ##  1   b0 0.1 0.21 0.4762 6.339e-01 -0.3116  0.5116
  ##  2   b1 1.3 0.32 4.0625 4.855e-05  0.6728  1.9272	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
