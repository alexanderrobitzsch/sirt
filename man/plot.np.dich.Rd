\name{plot.np.dich}
\alias{plot.np.dich}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Method for Object of Class \code{np.dich}
}
\description{
This function plots nonparametric item response
functions estimated with \code{dich.np}.
}
\usage{
\method{plot}{np.dich}(x, b, infit = NULL, outfit = NULL, 
    nsize = 100, askplot = TRUE, progress = TRUE, bands = FALSE, 
    plot.b = FALSE, shade = FALSE, shadecol = "burlywood1" , ...)		
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class \code{\link{np.dich}}
}
  \item{b}{
Estimated item difficulty (threshold)
}
  \item{infit}{
Infit (optional)
}
  \item{outfit}{
Outfit (optional)
}
  \item{nsize}{
XXX
}
  \item{askplot}{
Ask for new plot?
}
  \item{progress}{
Display progress?
}
  \item{bands}{
Draw confidence bands?
}
  \item{plot.b}{
Plot difficulty parameter?
}
  \item{shade}{
Shade curves?
}
  \item{shadecol}{
Shade color
}
  \item{\dots}{
Further arguments to be passed
}
}

\author{
Alexander Robitzsch
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For examples see \code{\link{np.dich}}.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Nonparametric item response theory}
\keyword{plot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
