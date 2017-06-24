
\name{lsdm}
\alias{lsdm}
\alias{summary.lsdm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Least Squares Distance Method of Cognitive Validation
}
\description{
This function estimates the least squares distance method
of cognitive validation (Dimitrov, 2007; Dimitrov & Atanasov, 2012) 
which assumes a multiplicative relationship of attribute response
probabilities to explain item response probabilities.
The function also estimates the classical linear logistic test 
model (LLTM; Fischer, 1973) which assumes a linear relationship
for item difficulties in the Rasch model.
}
\usage{
lsdm(data, Qmatrix, theta=qnorm(seq(5e-04,0.9995,len=100)), 
       quant.list=c(0.5,0.65,0.8), b=NULL, a=rep(1,nrow(Qmatrix)), 
       c=rep(0,nrow(Qmatrix)) )

\method{summary}{lsdm}(object,...)    
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
The discrete grid points where item response fuctions
are evaluated for doing the LSDM method.
}
  \item{quant.list}{
A vector of quantiles where attribute response functions are
evaluated.
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
  \item{object}{
Object of class \code{lsdm}
}
  \item{\dots}{
Further arguments to be passed
}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETAILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\details{
The least squares distance method (LSDM; Dimitrov 2007) is based on the
assumption that estimated item response functions \eqn{P(X_i = 1 | \theta)} 
can be decomposed in a multiplicative way (in the implemented 
conjunctive model):
\deqn{ P( X_i = 1  | \theta ) = \prod_{k=1}^K [ P( A_k = 1 | \theta ) ]^{q_{ik}} }
where \eqn{P( A_k = 1 | \theta )} are attribute response functions and
\eqn{q_{ik}} are entries of the Q-matrix. Note that the multiplicative form
can be rewritten by taking the logarithm
\deqn{ \log P( X_i = 1  | \theta ) = 
	\sum_{k=1}^K q_{ik} \log [ P( A_k = 1 | \theta ) ] }
Evaluation item and attribute response functions on a grid of \eqn{\theta} values and 
collecting these values in matrices \eqn{\bold{L}= \{ \log P( X_i = 1 ) | \theta ) \} }, 
\eqn{\bold{Q}= \{ q_{ik} \} } and
\eqn{\bold{X}= \{ \log P( A_k = 1  | \theta ) \} } leads to a least squares
problem of the form \eqn{ \bold{L} \approx \bold{Q} \bold{X}} with the
restriction of positive \bold{X} matrix entries. This least squares problem
is a linear inequality constrained model which is solved by making use of
the \pkg{ic.infer} package (Groemping, 2010).

After fitting the attribute response functions, empirical item-attribute
discriminations \eqn{w_{ik}} are calculated as the approximation of the following
equation
\deqn{ \log P( X_i = 1  | \theta ) = 
\sum_{k=1}^K w_{ik} q_{ik} \log [ P( A_k = 1 | \theta ) ] }
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
\item{data.fitted}{LSDM-fitted item reponse functions evaluated at \code{theta}}
\item{theta}{Grid of ability distributions at which 
functions are evaluated}
\item{item}{Item statistics (p value, \eqn{MAD}, ...)}
\item{data}{Estimated or fixed item reponse functions evaluated at \code{theta}}
\item{Qmatrix}{Used Q-matrix}
\item{lltm}{Model output of LLTM (\code{lm} values)}
\item{W}{Matrix with empirical item-attribute discriminations}
}
\references{
  DiBello, L. V., Roussos, L. A., & Stout, W. F. (2007). Review of
  cognitively diagnostic assessment and a summary of psychometric models.
  In C. R. Rao and S. Sinharay (Eds.), \emph{Handbook of Statistics},
  Vol. 26 (pp. 979-1030). Amsterdam: Elsevier.
  
  Dimitrov, D. M. (2007). Least squares distance method of cognitive validation 
  and analysis for binary items using their item response theory parameters. 
  \emph{Applied Psychological Measurement}, \bold{31}, 367-387.
   
  Dimitrov, D. M., & Atanasov, D. V. (2012). Conjunctive and disjunctive 
  extensions of the least squares distance model of cognitive diagnosis. 
  \emph{Educational and Psychological Measurement}, \bold{72}, 120-138.  
  
  Fischer, G. H. (1973). The linear logistic test model as an instrument 
  in educational research. \emph{Acta Psychologica}, \bold{37}, 359-374.
	
  Groemping, U. (2010). Inference with linear equality and inequality 
  constraints using \R: The package \pkg{ic.infer}.
  \emph{Journal of Statistical Software}, \bold{33(10)}, 1-31.  
  
Sonnleitner, P. (2008). Using the LLTM to evaluate an item-generating system 
for reading comprehension. \emph{Psychology Science}, \bold{50}, 345-362.  
}
\author{
Alexander Robitzsch
}
\note{
This function needs the \pkg{ic.infer} package.
}

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
# EXAMPLE 1: DATA FISCHER (see Dimitrov, 2007)
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
Qmatrix <- matrix( Qmatrix , nrow=29, byrow=TRUE )
colnames(Qmatrix) <- paste("A",1:8,sep="")
rownames(Qmatrix) <- paste("Item",1:29,sep="")

# Perform a LSDM analysis
lsdm.res <- sirt::lsdm( b = b, Qmatrix = Qmatrix )
summary(lsdm.res)
  ## Model Fit
  ## Model Fit LSDM   -  Mean MAD:  0.071     Median MAD:   0.07 
  ## Model Fit LLTM   -  Mean MAD:  0.079     Median MAD:  0.063    R^2= 0.615 
  ## ................................................................................ 
  ## Attribute Parameters
  ##    N.Items  b.2PL a.2PL  b.1PL eta.LLTM se.LLTM pval.LLTM
  ## A1      27 -2.101 1.615 -2.664   -1.168   0.404     0.009
  ## A2       8 -3.736 3.335 -5.491   -0.645   0.284     0.034
  ## A3      12 -5.491 0.360 -2.685   -0.013   0.284     0.963
  ## A4      22 -0.081 0.744 -0.059    1.495   0.350     0.000
  ## A5       7 -2.306 0.580 -1.622    0.243   0.301     0.428
  ## A6      10 -1.946 0.542 -1.306    0.447   0.243     0.080
  ## A7       5 -4.247 1.283 -4.799   -0.147   0.316     0.646
  ## A8       5 -2.670 0.663 -2.065    0.077   0.310     0.806
  ## [...]

#############################################################################
# EXAMPLE 2 DATA HENNING (see Dimitrov, 2007)
#############################################################################

# item difficulties
b <- c(-2.03,-1.29,-1.03,-1.58,0.59,-1.65,2.22,-1.46,2.58,-0.66)
# item slopes
a <- c(0.6,0.81,0.75,0.81,0.62,0.75,0.54,0.65,0.75,0.54)
# define Q-matrix
Qmatrix <- c(1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,1,1,0,0,
    0,0,0,1,0,0,1,0,0,1,0,0,0,1,0,0,0,0,1,1,1,0,1,0,0 )
Qmatrix <- matrix( Qmatrix , nrow=10, byrow=TRUE )
colnames(Qmatrix) <- paste("A",1:5,sep="")
rownames(Qmatrix) <- paste("Item",1:10,sep="")

# LSDM analysis
lsdm.res <- sirt::lsdm( b = b, a=a , Qmatrix = Qmatrix )
summary(lsdm.res)
  ## Model Fit LSDM   -  Mean MAD:  0.061     Median MAD:   0.06 
  ## Model Fit LLTM   -  Mean MAD:  0.069     Median MAD:  0.069    R^2= 0.902 
  ## ................................................................................ 
  ## Attribute Parameters
  ##    N.Items  b.2PL a.2PL  b.1PL eta.LLTM se.LLTM pval.LLTM
  ## A1       2 -2.727 0.786 -2.367   -1.592   0.478     0.021
  ## A2       5 -2.099 0.794 -1.834   -0.934   0.295     0.025
  ## A3       2 -0.763 0.401 -0.397    1.260   0.507     0.056
  ## A4       4 -1.459 0.638 -1.108   -0.738   0.309     0.062
  ## A5       2  2.410 0.509  1.564    2.673   0.451     0.002
  ## [...]
  ##
  ##  Discrimination Parameters
  ##  
  ##            A1    A2   A3    A4    A5
  ##  Item1  1.723    NA   NA    NA    NA
  ##  Item2     NA 1.615   NA    NA    NA
  ##  Item3     NA 0.650   NA 0.798    NA
  ##  Item4     NA 1.367   NA    NA    NA
  ##  Item5     NA 1.001 1.26    NA    NA
  ##  Item6     NA    NA   NA 0.866    NA
  ##  Item7     NA 0.697   NA    NA 0.891
  ##  Item8     NA    NA   NA 0.997    NA
  ##  Item9     NA    NA   NA 1.312 1.074
  ##  Item10 1.000    NA 0.74    NA    NA

\dontrun{
#############################################################################
# EXAMPLE 3: PISA reading (data.pisaRead)
#    using nonparametrically estimated item response functions
#############################################################################

data(data.pisaRead)
# response data
dat <- data.pisaRead$data
dat <- dat[ , substring( colnames(dat),1,1)=="R" ]
# define Q-matrix
pars <- data.pisaRead$item
Qmatrix <- data.frame(  "A0" = 1*(pars$ItemFormat=="MC" ) , 
                  "A1" = 1*(pars$ItemFormat=="CR" ) )

# start with estimating the 1PL in order to get person parameters
mod <- sirt::rasch.mml2( dat )
theta <- sirt::wle.rasch( dat=dat ,b = mod$item$b )$theta
# Nonparametric estimation of item response functions
mod2 <- sirt::np.dich( dat=dat , theta=theta , thetagrid = seq(-3,3,len=100) )

# LSDM analysis
lsdm.res <- sirt::lsdm( data=mod2$estimate , Qmatrix=Qmatrix , theta=mod2$thetagrid)
summary(lsdm.res)
  ## Model Fit
  ## Model Fit LSDM   -  Mean MAD:  0.215     Median MAD:  0.151 
  ## Model Fit LLTM   -  Mean MAD:  0.193     Median MAD:  0.119    R^2= 0.285 
  ## ................................................................................ 
  ## Attribute Parameter
  ##    N.Items  b.2PL a.2PL  b.1PL eta.LLTM se.LLTM pval.LLTM
  ## A0       5  1.326 0.705  1.289   -0.455   0.965     0.648
  ## A1       7 -1.271 1.073 -1.281   -1.585   0.816     0.081

#############################################################################
# EXAMPLE 4: Fraction subtraction dataset
#############################################################################

data( data.fraction1 , package="CDM")
data <- data.fraction1$data
q.matrix <- data.fraction1$q.matrix

#****
# Model 1: 2PL estimation 
mod1 <- sirt::rasch.mml2( data , est.a=1:nrow(q.matrix) )

# LSDM analysis
lsdm.res1 <- sirt::lsdm( b=mod1$item$b , a=mod1$item$a , Qmatrix=q.matrix )
summary(lsdm.res1)
  ##   
  ##   Model Fit LSDM   -  Mean MAD:  0.076     Median MAD:  0.055 
  ##   Model Fit LLTM   -  Mean MAD:  0.153     Median MAD:  0.155    R^2= 0.801 
  ##   ................................................................................ 
  ##   Attribute Parameter  
  ##       N.Items   b.2PL a.2PL  b.1PL eta.LLTM se.LLTM pval.LLTM
  ##   QT1      14  -0.741 2.991 -1.115   -0.815   0.217     0.004
  ##   QT2       8 -80.387 0.031 -4.806   -0.025   0.262     0.925
  ##   QT3      12  -2.461 0.711 -2.006   -0.362   0.268     0.207
  ##   QT4       9  -0.019 3.873 -0.100    1.465   0.345     0.002
  ##   QT5       3  -3.062 0.375 -1.481    0.254   0.280     0.387

#****
# Model 2: 1PL estimation  
mod2 <- sirt::rasch.mml2( data  )

# LSDM analysis
lsdm.res2 <- sirt::lsdm( b=mod1$item$b , Qmatrix=q.matrix )
summary(lsdm.res2)
  ##   
  ##   Model Fit LSDM   -  Mean MAD:  0.046     Median MAD:   0.03 
  ##   Model Fit LLTM   -  Mean MAD:  0.041     Median MAD:  0.042    R^2= 0.772

#############################################################################
# EXAMPLE 5: Dataset LLTM Sonnleitner Reading Comprehension (Sonnleitner, 2008)
#############################################################################

# item difficulties Table 7, p. 355 (Sonnleitner, 2008)
b <- c(-1.0189,1.6754,-1.0842,-.4457,-1.9419,-1.1513,2.0871,2.4874,-1.659,-1.197,-1.2437,
    2.1537,.3301,-.5181,-1.3024,-.8248,-.0278,1.3279,2.1454,-1.55,1.4277,.3301)
b <- b[-21] # remove Item 21

# Q-matrix Table 9 , p. 357 (Sonnleitner, 2008)
Qmatrix <- scan()    
   1 0 0 0 0 0 0 7 4 0 0 0   0 1 0 0 0 0 0 5 1 0 0 0   1 1 0 1 0 0 0 9 1 0 1 0   
   1 1 1 0 0 0 0 5 2 0 1 0   1 1 0 0 1 0 0 7 5 1 1 0   1 1 0 0 0 0 0 7 3 0 0 0   
   0 1 0 0 0 0 2 6 1 0 0 0   0 0 0 0 0 0 2 6 1 0 0 0   1 0 0 0 0 0 1 7 4 1 0 0   
   0 1 0 0 0 0 0 6 2 1 1 0   0 1 0 0 0 1 0 7 3 1 0 0   0 1 0 0 0 0 0 5 1 0 0 0   
   0 0 0 0 0 1 0 4 1 0 0 1   0 0 0 0 0 0 0 6 1 0 1 1   0 0 1 0 0 0 0 6 3 0 1 1   
   0 0 0 1 0 0 1 7 5 0 0 1   0 1 0 0 0 0 1 2 2 0 0 1   0 1 1 0 0 0 1 4 1 0 0 1   
   0 1 0 0 1 0 0 5 1 0 0 1   0 1 0 0 0 0 1 7 2 0 0 1   0 0 0 0 0 1 0 5 1 0 0 1
  
Qmatrix <- matrix( as.numeric(Qmatrix) , nrow=21 , ncol=12 , byrow=TRUE ) 
colnames(Qmatrix) <- scan( what="character" , nlines=1)
   pc ic ier inc iui igc ch nro ncro td a t

# divide Q-matrix entries by maximum in each column
Qmatrix <- round(Qmatrix / matrix(apply(Qmatrix,2,max),21,12,byrow=TRUE) ,3)
# LSDM analysis   
res <- sirt::lsdm( b=b , Qmatrix=Qmatrix )
summary(res)
  ##   
  ##   Model Fit LSDM   -  Mean MAD:  0.217     Median MAD:  0.178 
  ##   Model Fit LLTM   -  Mean MAD:  0.087     Median MAD:  0.062    R^2= 0.785 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Least Squares Distance Method (LSDM)}
\keyword{summary}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
