%% File Name: prmse.subscores.scales.Rd
%% File Version: 0.21
%% File Last Change: 2017-08-21 17:59:55

\name{prmse.subscores.scales}
\alias{prmse.subscores.scales}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Proportional Reduction of Mean Squared
Error (PRMSE) for Subscale Scores
}
\description{
This function estimates the proportional reduction of mean squared
error (PRMSE) according to
Haberman (Haberman 2008; Haberman, Sinharay & Puhan, 2008; see
Meijer et al. 2017 for an overview).
}
\usage{
prmse.subscores.scales(data, subscale)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
An \eqn{N \times I} data frame of item responses
}
  \item{subscale}{
Vector of labels corresponding to subscales
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Matrix with columns corresponding to subscales \cr
The symbol \code{X} denotes the subscale and \code{Z}
the whole scale (see also in the Examples section for the structure of
this matrix).
}
\references{
Haberman, S. J. (2008). When can subscores have value? 
\emph{Journal of Educational and Behavioral Statistics}, 
\bold{33}, 204-229.

Haberman, S., Sinharay, S., & Puhan, G. (2008).
Reporting subscores for institutions.
\emph{British Journal of Mathematical and Statistical Psychology},
\bold{62}, 79-95.

Meijer, R. R., Boeve, A. J., Tendeiro, J. N., Bosker, R. J., & Albers, C. J. (2017). 
The use of subscores in higher education: When is this useful?. 
\emph{Frontiers in Psychology | Educational Psychology, 8}.
}

\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See the \pkg{subscore} package for computing subscores and the PRMSE measures,
especially \code{subscore::CTTsub}.
}
%% \code{\link[subscore:CTTsub]{subscore::CTTsub}}

\examples{
#############################################################################
# EXAMPLE 1: PRMSE Reading data data.read
#############################################################################	

data( data.read )
p1 <- sirt::prmse.subscores.scales(data=data.read, 
         subscale = substring( colnames(data.read) , 1 ,1 ) )
print( p1 , digits= 3 )
  ##                 A       B       C
  ## N         328.000 328.000 328.000
  ## nX          4.000   4.000   4.000
  ## M.X         2.616   2.811   3.253
  ## Var.X       1.381   1.059   1.107
  ## SD.X        1.175   1.029   1.052
  ## alpha.X     0.545   0.381   0.640
  ## [...]
  ## nZ         12.000  12.000  12.000
  ## M.Z         8.680   8.680   8.680
  ## Var.Z       5.668   5.668   5.668
  ## SD.Z        2.381   2.381   2.381
  ## alpha.Z     0.677   0.677   0.677
  ## [...]
  ## cor.TX_Z    0.799   0.835   0.684
  ## rmse.X      0.585   0.500   0.505
  ## rmse.Z      0.522   0.350   0.614
  ## rmse.XZ     0.495   0.350   0.478
  ## prmse.X     0.545   0.381   0.640
  ## prmse.Z     0.638   0.697   0.468
  ## prmse.XZ    0.674   0.697   0.677
#-> Scales A and B do not have lower RMSE,
#   but for scale C the RMSE is smaller than the RMSE of a
#   prediction based on a whole scale.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Proportional reduction of mean squared error (PRMSE)}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
