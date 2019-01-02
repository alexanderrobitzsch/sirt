%% File Name: lsdm.Rd
%% File Version: 2.232


\name{lsdm}
\alias{lsdm}
\alias{summary.lsdm}
\alias{plot.lsdm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Least Squares Distance Method of Cognitive Validation
}
\description{
This function estimates the least squares distance method
of cognitive validation (Dimitrov, 2007; Dimitrov & Atanasov, 2012)
which assumes a multiplicative relationship of attribute response
probabilities to explain item response probabilities. The argument \code{distance}
allows the estimation of a squared loss function (\code{distance="L2"})
and an absolute value loss function (\code{distance="L1"}).

The function also estimates the classical linear logistic test
model (LLTM; Fischer, 1973) which assumes a linear relationship
for item difficulties in the Rasch model.
}

\usage{
lsdm(data, Qmatrix, theta=seq(-3,3,by=.5), wgt_theta=rep(1, length(theta)), distance="L2",
   quant.list=c(0.5,0.65,0.8), b=NULL, a=rep(1,nrow(Qmatrix)), c=rep(0,nrow(Qmatrix)) )

\method{summary}{lsdm}(object, file=NULL, digits=3, ...)

\method{plot}{lsdm}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
An \eqn{I \times L} matrix of dichotomous item responses.
The \code{data} consists of \eqn{I} item response functions
(parametrically or nonparametrically estimated) which are
evaluated at a discrete grid of \eqn{L} \code{theta} values
(person parameters)
and are specified in the argument \code{theta}.
 }
  \item{Qmatrix}{
An \eqn{I \times K} matrix where the allocation of items
to attributes is coded. Values of zero and one and all values between
zero and one are permitted. There must not be any items with only zero
Q-matrix entries in a row.
}
  \item{theta}{
The discrete grid points \eqn{\theta} where item response functions
are evaluated for doing the LSDM method.
}
\item{wgt_theta}{Optional vector for weights of discrete \eqn{\theta} points}
  \item{quant.list}{
A vector of quantiles where attribute response functions are
evaluated.
}
\item{distance}{Type of distance function for minimizing the discrepancy between
observed and expected item response functions. Options are \code{"L2"} which is the
squared distance (proposed in the original LSDM formulation in Dimitrov, 2007)
and the absolute value distance \code{"L1"} (see Details).
}
  \item{b}{
An optional vector of item difficulties. If it is specified,
then no \code{data} input is necessary.
}
  \item{a}{
An optional vector of item discriminations.
}
  \item{c}{
An optional vector of guessing parameters.
}
\item{object}{Object of class \code{lsdm}}
\item{file}{Optional file name for \code{summary} output}
\item{digits}{Number of digits aftert decimal in \code{summary}}
\item{\dots}{Further arguments to be passed}
\item{x}{Object of class \code{lsdm}}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETAILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\details{
The least squares distance method (LSDM; Dimitrov 2007) is based on the
assumption that estimated item response functions \eqn{P(X_i=1 | \theta)}
can be decomposed in a multiplicative way (in the implemented
conjunctive model):
\deqn{ P( X_i=1  | \theta ) \approx \prod_{k=1}^K [ P( A_k=1 | \theta ) ]^{q_{ik}} }
where \eqn{P( A_k=1 | \theta )} are attribute response functions and
\eqn{q_{ik}} are entries of the Q-matrix. Note that the multiplicative form
can be rewritten by taking the logarithm
\deqn{ \log P( X_i=1  | \theta ) \approx
    \sum_{k=1}^K q_{ik} \log [ P( A_k=1 | \theta ) ] }
The item and attribute response functions are evaluated on a grid of \eqn{\theta} values.
Using the definitions of matrices \eqn{\bold{L}=\{ \log P( X_i=1 ) | \theta ) \} },
\eqn{\bold{Q}=\{ q_{ik} \} } and
\eqn{\bold{X}=\{ \log P( A_k=1  | \theta ) \} }, the estimation problem can be formulated
as \eqn{ \bold{L} \approx \bold{Q} \bold{X}}. Two different loss functions for minimizing
the discrepancy between \eqn{ \bold{L}} and \eqn{\bold{Q} \bold{X}} are implemented.
First, the squared loss function computes the weighted difference
\eqn{|| \bold{L} - \bold{Q} \bold{X}||_2=\sum_i ( l_i - \sum_t q_{it} x_{it})^2}
 (\code{distance="L2"}) and has
been originally proposed by Dimitrov (2007). Second, the
absolute value loss function
\eqn{|| \bold{L} - \bold{Q} \bold{X}||_1=\sum_i | l_i - \sum_t q_{it} x_{it} |}
(\code{distance="L1"}) is more robust to outliers (i.e., items which
show misfit to the assumed multiplicative LSDM formulation).

After fitting the attribute response functions, empirical item-attribute
discriminations \eqn{w_{ik}} are calculated as the approximation of the following
equation
\deqn{ \log P( X_i=1  | \theta )=
\sum_{k=1}^K w_{ik} q_{ik} \log [ P( A_k=1 | \theta ) ] }
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\value{
A list with following entries
\item{mean.mad.lsdm0}{Mean of \eqn{MAD} statistics for LSDM}
\item{mean.mad.lltm}{Mean of \eqn{MAD} statistics for LLTM}
\item{attr.curves}{Estimated attribute response curves evaluated at \code{theta}}
\item{attr.pars}{Estimated attribute parameters for LSDM and LLTM}
\item{data.fitted}{LSDM-fitted item response functions evaluated at \code{theta}}
\item{theta}{Grid of ability distributions at which
functions are evaluated}
\item{item}{Item statistics (p value, \eqn{MAD}, ...)}
\item{data}{Estimated or fixed item response functions evaluated at \code{theta}}
\item{Qmatrix}{Used Q-matrix}
\item{lltm}{Model output of LLTM (\code{lm} values)}
\item{W}{Matrix with empirical item-attribute discriminations}
}


\references{

Al-Shamrani, A., & Dimitrov, D. M. (2016). Cognitive diagnostic analysis of reading
comprehension items: The case of English proficiency assessment in Saudi Arabia.
\emph{International Journal of School and Cognitive Psychology, 4}(3). 1000196.
http://dx.doi.org/10.4172/2469-9837.1000196

DiBello, L. V., Roussos, L. A., & Stout, W. F. (2007). Review of
cognitively diagnostic assessment and a summary of psychometric models.
In C. R. Rao and S. Sinharay (Eds.), \emph{Handbook of Statistics},
Vol. 26 (pp. 979-1030). Amsterdam: Elsevier.

Dimitrov, D. M. (2007). Least squares distance method of cognitive validation
and analysis for binary items using their item response theory parameters.
\emph{Applied Psychological Measurement, 31}, 367-387.
http://dx.doi.org/10.1177/0146621606295199

Dimitrov, D. M., & Atanasov, D. V. (2012). Conjunctive and disjunctive
extensions of the least squares distance model of cognitive diagnosis.
\emph{Educational and Psychological Measurement, 72}, 120-138.
http://dx.doi.org/10.1177/0013164411402324

Dimitrov, D. M., Gerganov, E. N., Greenberg, M., & Atanasov, D. V. (2008).
\emph{Analysis of cognitive attributes for mathematics items in the framework of
Rasch measurement}. AERA 2008, New York.

Fischer, G. H. (1973). The linear logistic test model as an instrument
in educational research. \emph{Acta Psychologica, 37}, 359-374.
http://dx.doi.org/10.1016/0001-6918(73)90003-6

Sonnleitner, P. (2008). Using the LLTM to evaluate an item-generating system
for reading comprehension. \emph{Psychology Science, 50}, 345-362.
}


%  Groemping, U. (2010). Inference with linear equality and inequality
%  constraints using \R: The package \pkg{ic.infer}.
%  \emph{Journal of Statistical Software, 33}(10), 1-31.

%\author{
%Alexander Robitzsch
%}


%\note{
%xx
%}


% This function needs the \pkg{ic.infer} package.

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Get a summary of the LSDM analysis with \code{\link{summary.lsdm}}.

See the \pkg{CDM} package for the estimation of related
cognitive diagnostic models (DiBello, Roussos & Stout, 2007).
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



\examples{
#############################################################################
# EXAMPLE 1: Dataset Fischer (see Dimitrov, 2007)
#############################################################################

# item difficulties
b <- c( 0.171,-1.626,-0.729,0.137,0.037,-0.787,-1.322,-0.216,1.802,
    0.476,1.19,-0.768,0.275,-0.846,0.213,0.306,0.796,0.089,
    0.398,-0.887,0.888,0.953,-1.496,0.905,-0.332,-0.435,0.346,
    -0.182,0.906)
# read Q-matrix
Qmatrix <- c( 1,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,
    1,0,1,1,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,1,1,0,0,1,0,1,0,1,0,0,0,
    1,0,1,0,1,1,0,0,1,0,1,1,0,1,0,0,1,0,0,1,0,1,0,0,1,0,1,1,1,0,0,0,
    1,0,0,1,0,0,1,0,1,0,0,1,0,0,1,0,1,0,1,0,0,0,1,0,1,1,0,1,0,1,1,0,
    1,0,1,1,0,0,1,0,1,0,0,1,0,0,0,1,1,0,1,1,0,0,0,1,1,0,0,1,0,0,0,1,
    0,1,0,0,0,1,0,1,1,1,0,1,0,1,0,1,1,0,0,1,0,1,0,0,1,1,0,0,1,0,0,0,
    1,0,0,1,1,0,0,0,1,1,0,1,0,0,0,0,1,0,1,1,0,0,0,0,1,0,1,1,0,1,0,0,
    1,1,0,1,0,0,0,0,1,0,1,1,1,1,0,0 )
Qmatrix <- matrix( Qmatrix, nrow=29, byrow=TRUE )
colnames(Qmatrix) <- paste("A",1:8,sep="")
rownames(Qmatrix) <- paste("Item",1:29,sep="")

#* Model 1: perform a LSDM analysis with defaults
mod1 <- sirt::lsdm( b=b, Qmatrix=Qmatrix )
summary(mod1)
plot(mod1)

#* Model 2: different theta values and weights
theta <- seq(-4,4,len=31)
wgt_theta <- stats::dnorm(theta)
mod2 <- sirt::lsdm( b=b, Qmatrix=Qmatrix, theta=theta, wgt_theta=wgt_theta )
summary(mod2)

#* Model 3: absolute value distance function
mod3 <- sirt::lsdm( b=b, Qmatrix=Qmatrix, distance="L1" )
summary(mod3)

#############################################################################
# EXAMPLE 2: Dataset Henning (see Dimitrov, 2007)
#############################################################################

# item difficulties
b <- c(-2.03,-1.29,-1.03,-1.58,0.59,-1.65,2.22,-1.46,2.58,-0.66)
# item slopes
a <- c(0.6,0.81,0.75,0.81,0.62,0.75,0.54,0.65,0.75,0.54)
# define Q-matrix
Qmatrix <- c(1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,1,1,0,0,
    0,0,0,1,0,0,1,0,0,1,0,0,0,1,0,0,0,0,1,1,1,0,1,0,0 )
Qmatrix <- matrix( Qmatrix, nrow=10, byrow=TRUE )
colnames(Qmatrix) <- paste("A",1:5,sep="")
rownames(Qmatrix) <- paste("Item",1:10,sep="")

# LSDM analysis
mod <- sirt::lsdm( b=b, a=a, Qmatrix=Qmatrix )
summary(mod)

\dontrun{
#############################################################################
# EXAMPLE 3: PISA reading (data.pisaRead)
#    using nonparametrically estimated item response functions
#############################################################################

data(data.pisaRead)
# response data
dat <- data.pisaRead$data
dat <- dat[, substring( colnames(dat),1,1)=="R" ]
# define Q-matrix
pars <- data.pisaRead$item
Qmatrix <- data.frame(  "A0"=1*(pars$ItemFormat=="MC" ),
                  "A1"=1*(pars$ItemFormat=="CR" ) )

# start with estimating the 1PL in order to get person parameters
mod <- sirt::rasch.mml2( dat )
theta <- sirt::wle.rasch( dat=dat,b=mod$item$b )$theta
# Nonparametric estimation of item response functions
mod2 <- sirt::np.dich( dat=dat, theta=theta, thetagrid=seq(-3,3,len=100) )

# LSDM analysis
lmod <- sirt::lsdm( data=mod2$estimate, Qmatrix=Qmatrix, theta=mod2$thetagrid)
summary(lmod)
plot(lmod)

#############################################################################
# EXAMPLE 4: Fraction subtraction dataset
#############################################################################

data( data.fraction1, package="CDM")
data <- data.fraction1$data
q.matrix <- data.fraction1$q.matrix

#****
# Model 1: 2PL estimation
mod1 <- sirt::rasch.mml2( data, est.a=1:nrow(q.matrix) )

# LSDM analysis
lmod1 <- sirt::lsdm( b=mod1$item$b, a=mod1$item$a, Qmatrix=q.matrix )
summary(lmod1)

#****
# Model 2: 1PL estimation
mod2 <- sirt::rasch.mml2(data)

# LSDM analysis
lmod2 <- sirt::lsdm( b=mod1$item$b, Qmatrix=q.matrix )
summary(lmod2)

#############################################################################
# EXAMPLE 5: Dataset LLTM Sonnleitner Reading Comprehension (Sonnleitner, 2008)
#############################################################################

# item difficulties Table 7, p. 355 (Sonnleitner, 2008)
b <- c(-1.0189,1.6754,-1.0842,-.4457,-1.9419,-1.1513,2.0871,2.4874,-1.659,-1.197,-1.2437,
    2.1537,.3301,-.5181,-1.3024,-.8248,-.0278,1.3279,2.1454,-1.55,1.4277,.3301)
b <- b[-21] # remove Item 21

# Q-matrix Table 9, p. 357 (Sonnleitner, 2008)
Qmatrix <- scan()
   1 0 0 0 0 0 0 7 4 0 0 0   0 1 0 0 0 0 0 5 1 0 0 0   1 1 0 1 0 0 0 9 1 0 1 0
   1 1 1 0 0 0 0 5 2 0 1 0   1 1 0 0 1 0 0 7 5 1 1 0   1 1 0 0 0 0 0 7 3 0 0 0
   0 1 0 0 0 0 2 6 1 0 0 0   0 0 0 0 0 0 2 6 1 0 0 0   1 0 0 0 0 0 1 7 4 1 0 0
   0 1 0 0 0 0 0 6 2 1 1 0   0 1 0 0 0 1 0 7 3 1 0 0   0 1 0 0 0 0 0 5 1 0 0 0
   0 0 0 0 0 1 0 4 1 0 0 1   0 0 0 0 0 0 0 6 1 0 1 1   0 0 1 0 0 0 0 6 3 0 1 1
   0 0 0 1 0 0 1 7 5 0 0 1   0 1 0 0 0 0 1 2 2 0 0 1   0 1 1 0 0 0 1 4 1 0 0 1
   0 1 0 0 1 0 0 5 1 0 0 1   0 1 0 0 0 0 1 7 2 0 0 1   0 0 0 0 0 1 0 5 1 0 0 1

Qmatrix <- matrix( as.numeric(Qmatrix), nrow=21, ncol=12, byrow=TRUE )
colnames(Qmatrix) <- scan( what="character", nlines=1)
   pc ic ier inc iui igc ch nro ncro td a t

# divide Q-matrix entries by maximum in each column
Qmatrix <- round(Qmatrix / matrix(apply(Qmatrix,2,max),21,12,byrow=TRUE),3)
# LSDM analysis
mod <- sirt::lsdm( b=b, Qmatrix=Qmatrix )
summary(mod)

#############################################################################
# EXAMPLE 6: Dataset Dimitrov et al. (2008)
#############################################################################

Qmatrix <- scan()
1 0 0 0 1 1 0 1 1 0 0 0 1 0 0 0 1 1 0 0 1 1 0 1 0 0 1 1 1 0 1 0 0 0 1 0

Qmatrix <- matrix(Qmatrix, ncol=4, byrow=TRUE)
colnames(Qmatrix) <- paste0("A",1:4)
rownames(Qmatrix) <- paste0("I",1:9)

b <- scan()
0.068 1.095 -0.641 -1.129 -0.061 1.218 1.244 -0.648 -1.146

# estimate model
mod <- sirt::lsdm( b=b, Qmatrix=Qmatrix )
summary(mod)
plot(mod)

#############################################################################
# EXAMPLE 7: Dataset Al-Shamrani & Dimitrov et al. (2017)
#############################################################################

I <- 39  # number of items

Qmatrix <- scan()
0 0 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0
0 0 0 0 0 1 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0
0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0
0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0
0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0
0 1 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 1 0 0 0
0 0 1 0 0 0 0 0 0 0 0 0 0 1 0

Qmatrix <- matrix(Qmatrix, nrow=I, byrow=TRUE)
colnames(Qmatrix) <- paste0("A",1:7)
rownames(Qmatrix) <- paste0("I",1:I)

pars <- scan()
1.952 0.9833 0.1816 1.1053 0.9631 0.1653 1.3904 1.3208 0.2545 0.7391 1.9367 0.2083 2.0833
1.8627 0.1873 1.4139 1.0107 0.2454 0.8274 0.9913 0.2137 1.0338 -0.0068 0.2368 2.4803
0.7939 0.1997 1.4867 1.1705 0.2541 1.4482 1.4176 0.2889 1.0789 0.8062 0.269 1.6258 1.1739
0.1723 1.5995 1.0936 0.2054 1.1814 1.0909 0.2623 2.0389 1.5023 0.2466 1.3636 1.1485 0.2059
1.8468 1.2755 0.192 1.9461 1.4947 0.2001 1.194 0.0889 0.2275 1.2114 0.8925 0.2367 2.0912
0.5961 0.2036 2.5769 1.3014 0.186 1.4554 1.2529 0.2423 1.4919 0.4763 0.2482 2.6787 1.7069
0.1796 1.5611 1.3991 0.2312 1.4353 0.678 0.1851 0.9127 1.3523 0.2525 0.6886 -0.3652 0.207
0.7039 -0.2494 0.2315 1.3683 0.8953 0.2326 1.4992 0.1025 0.2403 1.0727 0.2591 0.2152
1.3854 1.3802 0.2448 0.7748 0.4304 0.184 1.0218 1.8964 0.1949 1.5773 1.8934 0.2231 0.8631
1.4145 0.2132

pars <- matrix(pars, nrow=I, byrow=TRUE)
colnames(pars) <- c("a","b","c")
rownames(pars) <- paste0("I",1:I)
pars <- as.data.frame(pars)

#* Model 1: fit LSDM to 3PL curves (as in Al-Shamrani)
mod1 <- sirt::lsdm(b=pars$b, a=pars$a, c=pars$c, Qmatrix=Qmatrix)
summary(mod1)
plot(mod1)

#* Model 2: fit LSDM to 2PL curves
mod2 <- sirt::lsdm(b=pars$b, a=pars$a, Qmatrix=Qmatrix)
summary(mod2)
plot(mod2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Least Squares Distance Method (LSDM)}
%% \keyword{summary}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
