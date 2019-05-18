%% File Name: data.si.Rd
%% File Version: 0.395

\name{data.sirt}

\alias{data.sirt}
\alias{data.si01}
\alias{data.si02}
\alias{data.si03}
\alias{data.si04}
\alias{data.si05}
\alias{data.si06}
\alias{data.si07}
\alias{data.si08}


\docType{data}


\title{
Some Example Datasets for the \code{sirt} Package
}


\description{
Some example datasets for the \code{sirt} package.
}


\usage{
data(data.si01)
data(data.si02)
data(data.si03)
data(data.si04)
data(data.si05)
data(data.si06)
data(data.si07)
data(data.si08)
}


\format{

\itemize{
\item The format of the dataset \code{data.si01} is:

\code{'data.frame':   1857 obs. of  3 variables:} \cr
\code{ $ idgroup: int  1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{ $ item1  : int  NA NA NA NA NA NA NA NA NA NA ...} \cr
\code{ $ item2  : int  4 4 4 4 4 4 4 2 4 4 ...} \cr

\item The dataset \code{data.si02} is the Stouffer-Toby-dataset published
in Lindsay, Clogg and Grego (1991; Table 1, p.97, Cross-classification A):

\code{List of 2} \cr
\code{ $ data   : num [1:16, 1:4] 1 0 1 0 1 0 1 0 1 0 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 2} \cr
\code{  .. ..$ : NULL} \cr
\code{  .. ..$ : chr [1:4] "I1" "I2" "I3" "I4"} \cr
\code{ $ weights: num [1:16] 42 1 6 2 6 1 7 2 23 4 ...} \cr

\item The format of the dataset \code{data.si03} (containing item
parameters of two studies) is:

\code{'data.frame':   27 obs. of  3 variables:} \cr
\code{ $ item    : Factor w/ 27 levels "M1","M10","M11",..: 1 12 21 22 ...} \cr
\code{ $ b_study1: num  0.297 1.163 0.151 -0.855 -1.653 ...} \cr
\code{ $ b_study2: num  0.72 1.118 0.351 -0.861 -1.593 ...} \cr

\item The dataset \code{data.si04} is adapted from Bartolucci, Montanari
and Pandolfi (2012; Table 4, Table 7). The data contains 4999 persons,
79 items on 5 dimensions. See \code{\link{rasch.mirtlc}} for using the
data in an analysis.

\code{List of 3} \cr
\code{ $ data        : num [1:4999, 1:79] 0 1 1 0 1 1 0 0 1 1 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 2} \cr
\code{  .. ..$ : NULL} \cr
\code{  .. ..$ : chr [1:79] "A01" "A02" "A03" "A04" ...} \cr
\code{ $ itempars    :'data.frame':   79 obs. of  4 variables:} \cr
\code{  ..$ item      : Factor w/ 79 levels "A01","A02","A03",..: 1 2 3 4 5 6 7 8 9 10 ...} \cr
\code{  ..$ dim       : num [1:79] 1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{  ..$ gamma     : num [1:79] 1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{  ..$ gamma.beta: num [1:79] -0.189 0.25 0.758 1.695 1.022 ...} \cr
\code{ $ distribution: num [1:9, 1:7] 1 2 3 4 5 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 2} \cr
\code{  .. ..$ : NULL} \cr
\code{  .. ..$ : chr [1:7] "class" "A" "B" "C" ...} \cr

\item The dataset \code{data.si05} contains double ratings of two
exchangeable raters for three items which are in \code{Ex1}, \code{Ex2}
and \code{Ex3}, respectively.

\code{List of 3} \cr
\code{ $ Ex1:'data.frame':    199 obs. of  2 variables:} \cr
\code{  ..$ C7040: num [1:199] NA 1 0 1 1 0 0 0 1 0 ...} \cr
\code{  ..$ C7041: num [1:199] 1 1 0 0 0 0 0 0 1 0 ...} \cr
\code{ $ Ex2:'data.frame':    2000 obs. of  2 variables:} \cr
\code{  ..$ rater1: num [1:2000] 2 0 3 1 2 2 0 0 0 0 ...} \cr
\code{  ..$ rater2: num [1:2000] 4 1 3 2 1 0 0 0 0 2 ...} \cr
\code{ $ Ex3:'data.frame':    2000 obs. of  2 variables:} \cr
\code{  ..$ rater1: num [1:2000] 5 1 6 2 3 3 0 0 0 0 ...} \cr
\code{  ..$ rater2: num [1:2000] 7 2 6 3 2 1 0 1 0 3 ...} \cr

\item The dataset \code{data.si06} contains multiple choice item
responses. The correct alternative is denoted as 0, distractors
are indicated by the codes 1, 2 or 3.

\code{'data.frame':   4441 obs. of  14 variables:} \cr
\code{ $ WV01: num  0 0 0 0 0 0 0 0 0 3 ...} \cr
\code{ $ WV02: num  0 0 0 3 0 0 0 0 0 1 ...} \cr
\code{ $ WV03: num  0 1 0 0 0 0 0 0 0 0 ...} \cr
\code{ $ WV04: num  0 0 0 0 0 0 0 0 0 1 ...} \cr
\code{ $ WV05: num  3 1 1 1 0 0 1 1 0 2 ...} \cr
\code{ $ WV06: num  0 1 3 0 0 0 2 0 0 1 ...} \cr
\code{ $ WV07: num  0 0 0 0 0 0 0 0 0 0 ...} \cr
\code{ $ WV08: num  0 1 1 0 0 0 0 0 0 0 ...} \cr
\code{ $ WV09: num  0 0 0 0 0 0 0 0 0 2 ...} \cr
\code{ $ WV10: num  1 1 3 0 0 2 0 0 0 0 ...} \cr
\code{ $ WV11: num  0 0 0 0 0 0 0 0 0 0 ...} \cr
\code{ $ WV12: num  0 0 0 2 0 0 2 0 0 0 ...} \cr
\code{ $ WV13: num  3 1 1 3 0 0 3 0 0 0 ...} \cr
\code{ $ WV14: num  3 1 2 3 0 3 1 3 3 0 ...} \cr

\item The dataset \code{data.si07} contains parameters of the empirical illustration
of DeCarlo (XXXX). The simulation function \code{sim_fun} can be used for
simulating data from the IRSDT model (see DeCarlo, XXXX)

\code{List of 3} \cr
\code{ $ pars   :'data.frame':        16 obs. of  3 variables:} \cr
\code{  ..$ item: Factor w/ 16 levels "I01","I02","I03",..: 1 2 3 4 5 6 7 8 9 10 ...} \cr
\code{  ..$ b   : num [1:16] -1.1 -0.18 1.44 1.78 -1.19 0.45 -1.12 0.33 0.82 -0.43 ...} \cr
\code{  ..$ d   : num [1:16] 2.69 4.6 6.1 3.11 3.2 ...} \cr
\code{ $ trait  :'data.frame':        20 obs. of  2 variables:} \cr
\code{  ..$ x   : num [1:20] 0.025 0.075 0.125 0.175 0.225 0.275 0.325 0.375 0.425 0.475 ...} \cr
\code{  ..$ prob: num [1:20] 0.0238 0.1267 0.105 0.0594 0.0548 ...} \cr
\code{ $ sim_fun:function (lambda, b, d, items)  } \cr

\item The dataset \code{data.si08} contains 5 items with respect to knowledge
about lung cancer and the kind of information acquisition (Goodman, 1970;
see also Rasch, Kubinger & Yanagida, 2011).
\code{L1}: reading newspapers, \code{L2}: listening radio,
\code{L3}: reading books and magazines,
\code{L4}: attending talks, \code{L5}: knowledge about lung cancer

\code{'data.frame':   32 obs. of  6 variables:} \cr
\code{ $ L1 : num  1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{ $ L2 : num  1 1 1 1 1 1 1 1 0 0 ...} \cr
\code{ $ L3 : num  1 1 1 1 0 0 0 0 1 1 ...} \cr
\code{ $ L4 : num  1 1 0 0 1 1 0 0 1 1 ...} \cr
\code{ $ L5 : num  1 0 1 0 1 0 1 0 1 0 ...} \cr
\code{ $ wgt: num  23 8 102 67 8 4 35 59 27 18 ...} \cr
}
}





%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}


\seealso{
Some free datasets can be obtained from \cr
Psychological questionnaires: http://personality-testing.info/_rawdata/ \cr
PISA 2012: http://pisa2012.acer.edu.au/downloads.php \cr
PIAAC: http://www.oecd.org/site/piaac/publicdataandanalysis.htm \cr
TIMSS 2011: http://timssandpirls.bc.edu/timss2011/international-database.html \cr
ALLBUS: http://www.gesis.org/allbus/allbus-home/
}

%\source{
%  ~~ reference to a publication or URL from which the data were obtained ~~
%}



\references{
Bartolucci, F., Montanari, G. E., & Pandolfi, S. (2012).
Dimensionality of the latent structure and item selection via latent
class multidimensional IRT models. \emph{Psychometrika, 77}, 782-802.
\doi{10.1007/s11336-012-9278-0}

DeCarlo, L. T. (XXXX). An item response model for true-false exams based on signal
detection theory. \emph{Applied Psychological Measurement (online first)}.
\doi{10.1177/0146621619843823}

Goodman, L. A. (1970). The multivariate analysis of qualitative data: Interactions
among multiple classifications.
\emph{Journal of the American Statistical Association, 65}(329), 226-256.
\doi{10.1080/01621459.1970.10481076}

Lindsay, B., Clogg, C. C., & Grego, J. (1991).
Semiparametric estimation in the Rasch model and related exponential response
models, including a simple latent class model for item analysis.
\emph{Journal of the American Statistical Association, 86}, 96-107.
\doi{10.1080/01621459.1991.10475008}

Rasch, D., Kubinger, K. D., & Yanagida, T. (2011).
\emph{Statistics in psychology using R and SPSS}. New York: Wiley.
\doi{10.1002/9781119979630}
}


\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Nested logit model multiple choice dataset data.si06
#############################################################################

data(data.si06, package="sirt")
dat <- data.si06

#** estimate 2PL nested logit model
library(mirt)
mod1 <- mirt::mirt( dat, model=1, itemtype="2PLNRM", key=rep(0,ncol(dat) ),
            verbose=TRUE  )
summary(mod1)
cmod1 <- sirt::mirt.wrapper.coef(mod1)$coef
cmod1[,-1] <- round( cmod1[,-1], 3)

#** normalize item parameters according Suh and Bolt (2010)
cmod2 <- cmod1

# slope parameters
ind <-  grep("ak",colnames(cmod2))
h1 <- cmod2[,ind ]
cmod2[,ind] <- t( apply( h1, 1, FUN=function(ll){ ll - mean(ll) } ) )
# item intercepts
ind <-  paste0( "d", 0:9 )
ind <- which( colnames(cmod2) \%in\% ind )
h1 <- cmod2[,ind ]
cmod2[,ind] <- t( apply( h1, 1, FUN=function(ll){ ll - mean(ll) } ) )
cmod2[,-1] <- round( cmod2[,-1], 3)

#############################################################################
# EXAMPLE 2: Item response modle based on signal detection theory (IRSDT model)
#############################################################################

data(data.si07, package="sirt")
data <- data.si07

#-- simulate data
set.seed(98)
N <- 2000 # define sample size
# generate membership scores
lambda <- sample(size=N, x=data$trait$x, prob=data$trait$prob, replace=TRUE)
b <- data$pars$b
d <- data$pars$d
items <- data$pars$item
dat <- data$sim_fun(lambda=lambda, b=b, d=d, items=items)

#- estimate IRSDT model as a grade of membership model with two classes
problevels <- seq( 0.025, 0.975, length=20 )
mod1 <- sirt::gom.em( dat, K=2, problevels=problevels )
summary(mod1)
}
}
%% \keyword{datasets}
