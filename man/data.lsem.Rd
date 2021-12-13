%% File Name: data.lsem.Rd
%% File Version: 0.094

\name{data.lsem}
\alias{data.lsem01}
\alias{data.lsem02}
\alias{data.lsem03}

\docType{data}

\title{
Datasets for Local Structural Equation Models / Moderated Factor Analysis
}


\description{
Datasets for local structural equation models or moderated factor analysis.
}


\usage{
data(data.lsem01)
data(data.lsem02)
data(data.lsem03)
}



\format{

\itemize{
\item The dataset \code{data.lsem01} has the following structure

\code{'data.frame':   989 obs. of  6 variables:} \cr
\code{ $ age: num  4 4 4 4 4 4 4 4 4 4 ...} \cr
\code{ $ v1 : num  1.83 2.38 1.85 4.53 -0.04 4.35 2.38 1.83 4.81 2.82 ...} \cr
\code{ $ v2 : num  6.06 9.08 7.41 8.24 6.18 7.4 6.54 4.28 6.43 7.6 ...} \cr
\code{ $ v3 : num  1.42 3.05 6.42 -1.05 -1.79 4.06 -0.17 -2.64 0.84 6.42 ...} \cr
\code{ $ v4 : num  3.84 4.24 3.24 3.36 2.31 6.07 4 5.93 4.4 3.49 ...} \cr
\code{ $ v5 : num  7.84 7.51 6.62 8.02 7.12 7.99 7.25 7.62 7.66 7.03 ...} \cr

\item The dataset \code{data.lsem02} is a slightly perturbed dataset of the
Woodcock-Johnson  III (WJ-III) Tests of Cognitive Abilities used in Hildebrandt et al.
(2016) and has the following structure

\code{'data.frame':   1129 obs. of  8 variables:} \cr
\code{ $ age : int  4 4 4 4 4 4 4 4 4 4 ...} \cr
\code{ $ gcw : num  -3.53 -3.73 -3.77 -3.84 -4.26 -4.6 -3.66 -4.31 -4.46 -3.64 ...} \cr
\code{ $ gvw : num  -1.98 -1.35 -1.66 -3.24 -1.17 -2.78 -2.97 -3.88 -3.22 -0.68 ...} \cr
\code{ $ gfw : num  -2.49 -2.41 -4.48 -4.17 -4.43 -5.06 -3.94 -3.66 -3.7 -2.74 ...} \cr
\code{ $ gsw : num  -4.85 -5.05 -5.66 -4.3 -5.23 -5.63 -4.91 -5.75 -6.29 -5.47 ...} \cr
\code{ $ gsmw: num  -2.99 -1.13 -4.21 -3.59 -3.79 -4.77 -2.98 -4.48 -2.99 -3.83 ...} \cr
\code{ $ glrw: num  -2.49 -2.91 -3.45 -2.91 -3.31 -3.78 -3.5 -3.96 -2.97 -3.14 ...} \cr
\code{ $ gaw : num  -3.22 -3.77 -3.54 -3.6 -3.22 -3.5 -1.27 -2.08 -2.23 -3.25 ...} \cr

\item The dataset \code{data.lsem03} is a synthetic dataset of the SON-R application
used in Hueluer et al. (2011) has the following structure

\code{'data.frame':   1027 obs. of  10 variables:} \cr
\code{ $ id       : num  10001 10002 10003 10004 10005 ...} \cr
\code{ $ female   : int  0 0 0 0 0 0 0 0 0 0 ...} \cr
\code{ $ age      : num  2.62 2.65 2.66 2.67 2.68 2.68 2.68 2.69 2.71 2.71 ...} \cr
\code{ $ age_group: int  1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{ $ p1       : num  -1.98 -1.98 -1.67 -2.29 -1.67 -1.98 -2.29 -1.98 -2.6 -1.67 ...} \cr
\code{ $ p2       : num  -1.51 -1.51 -0.55 -1.84 -1.51 -1.84 -2.16 -1.84 -2.48 -1.84 ...} \cr
\code{ $ p3       : num  -1.4 -2.31 -1.1 -2 -1.4 -1.7 -2.31 -1.4 -2.31 -0.79 ...} \cr
\code{ $ r1       : num  -1.46 -1.14 -0.49 -2.11 -1.46 -1.46 -2.11 -1.46 -2.75 -1.78 ...} \cr
\code{ $ r2       : num  -2.67 -1.74 0.74 -1.74 -0.81 -1.43 -2.05 -1.43 -1.74 -1.12 ...} \cr
\code{ $ r3       : num  -1.64 -1.64 -1.64 -0.9 -1.27 -3.11 -2.74 -1.64 -2.37 -1.27 ...} \cr

The subtests  Mosaics (\code{p1}), Puzzles (\code{p1}), and Patterns (\code{p3})
constitute the performance subscale;
the subtests Categories (\code{r1}), Analogies (\code{r2}), and
Situations (\code{r3}) constitute the reasoning subscale.
}
}


%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{
Hildebrandt, A., Luedtke, O., Robitzsch, A., Sommer, C., &
Wilhelm, O. (2016). Exploring factor model parameters across continuous variables
with local structural equation models.
\emph{Multivariate Behavioral Research, 51}(2-3), 257-278.
\doi{10.1080/00273171.2016.1142856}

Hueluer, G., Wilhelm, O., & Robitzsch, A. (2011). Intelligence differentiation in
early childhood. \emph{Journal of Individual Differences, 32}(3), 170-179.
\doi{10.1027/1614-0001/a000049}
}

%\examples{
%data(data.lsem01)
%}

%% \keyword{datasets}
