%% File Name: linking_2groups.Rd
%% File Version: 0.054

\name{linking_2groups}
\alias{linking_2groups}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linking Two Groups
}
\description{
This function performs (robust) Haebara linking or Stocking-Lord for two groups
in the 2PL model. There is also an option for simultaneous estimation in which
joint item parameters should also be estimated.
}
\usage{
linking_2groups(pars, method, type="asymm", pow=2, eps=0.001, simultan=FALSE,
     Theta=seq(-6, 6, len=101), wgt=NULL, par_init=NULL, optimizer="nlminb",
     control_optimizer=list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pars}{
Matrix containing item parameters that must have columns with
names \code{a1}, \code{b1}, \code{a2}, and \code{b2}.
}
  \item{method}{
Linking method. Can be \code{SL} (Stocking-Lord) or \code{Hae} (Haebara).
}
  \item{type}{
Type of linking function. Can be \code{"asymm"} or \code{"symm"} for
asymmetric or symmetric linking, respectively.
}
  \item{pow}{
Power used in the loss function. The values should be larger than or
equal than zero. The default values is 2.
}
  \item{eps}{
Tuning parameter in numerical differentiation
}
  \item{simultan}{
Logical indicating whether joint item parameters should be simultaneously
estimated
}
  \item{Theta}{
Theta grid
}
  \item{wgt}{
Weights for theta grid
}
  \item{par_init}{
Initial vector for parameters
}
  \item{optimizer}{
Chosen optimizer. Can be \code{"nlminb"} or \code{"optim"}.
}
  \item{control_optimizer}{
Control arguments for optimizer
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
List with following entries
\item{par}{Estimated linking parameters}
\item{$\ldots$}{More values}
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: data.ex16 from TAM package
#############################################################################

#- extract dataset
data(data.ex16, package="TAM")
dat <- data.ex16
items <- colnames(dat)[-c(1,2)]

# fit grade 1
rdat1 <- TAM::tam_remove_missings( dat[ dat$grade==1, ], items=items )
mod1 <- sirt::rasch.mml2( dat=rdat1$resp[, rdat1$items], est.a=1:length(rdat1$items) )

# fit grade 2
rdat2 <- TAM::tam_remove_missings( dat[ dat$grade==2, ], items=items )
mod2 <- sirt::rasch.mml2( dat=rdat2$resp[, rdat2$items], est.a=1:length(rdat2$items) )

# merge item parameters of separate scalings
item1 <- mod1$item[, c("item", "a","b") ]
colnames(item1)[2:3] <- c("a1","b1")
item2 <- mod2$item[, c("item", "a","b") ]
colnames(item2)[2:3] <- c("a2","b2")
pars <- merge(x=item1, y=item2, by="item")

# choose theta grid
Theta <- seq(-6,6, len=101)
wgt <- sirt::sirt_dnorm_discrete(x=Theta, mean=0, sd=2)

#- robust Haebara linking with powers 0.5 or 0 and different types
res1 <- sirt::linking_2groups(pars=pars, method="Hae", type="asymm",
              Theta=Theta, wgt=wgt, pow=0.5, eps=0.001 )
res2 <- sirt::linking_2groups(pars=pars, method="Hae", type="symm",
              Theta=Theta, wgt=wgt, pow=0, eps=0.01 )
res1$par  # linking parameter estimate
res2$par
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
