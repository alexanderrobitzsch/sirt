%% File Name: data.bs.Rd
%% File Version: 0.08


\name{data.bs}
\alias{data.bs}
\alias{data.bs07a}

\docType{data}


\title{
Datasets from Borg and Staufenbiel (2007)
}

\description{
Datasets of the book of Borg and Staufenbiel (2007)
\emph{Lehrbuch Theorien and Methoden der Skalierung}.
}


\usage{
data(data.bs07a)
}


\format{
\itemize{
\item The dataset \code{data.bs07a} contains the data
\emph{Gefechtsangst} (p. 130) and contains 8 of the original 9 items.
The items are symptoms of anxiety in engagement. \cr
\code{GF1}: starkes Herzklopfen, \code{GF2}: flaues Gefuehl in der
Magengegend, \code{GF3}: Schwaechegefuehl, \code{GF4}: Uebelkeitsgefuehl,
\code{GF5}: Erbrechen, \code{GF6}: Schuettelfrost,
\code{GF7}: in die Hose urinieren/einkoten, \code{GF9}: Gefuehl der
Gelaehmtheit

The format is

\code{'data.frame':   100 obs. of  9 variables:} \cr
\code{ $ idpatt: int  44 29 1 3 28 50 50 36 37 25 ...} \cr
\code{ $ GF1   : int  1 1 1 1 1 0 0 1 1 1 ...} \cr
\code{ $ GF2   : int  0 1 1 1 1 0 0 1 1 1 ...} \cr
\code{ $ GF3   : int  0 0 1 1 0 0 0 0 0 1 ...} \cr
\code{ $ GF4   : int  0 0 1 1 0 0 0 1 0 1 ...} \cr
\code{ $ GF5   : int  0 0 1 1 0 0 0 0 0 0 ...} \cr
\code{ $ GF6   : int  1 1 1 1 1 0 0 0 0 0 ...} \cr
\code{ $ GF7   : num  0 0 1 1 0 0 0 0 0 0 ...} \cr
\code{ $ GF9   : int  0 0 1 1 1 0 0 0 0 0 ...} \cr

\item \emph{MORE DATASETS}

  }
}


%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}


\references{
Borg, I., & Staufenbiel, T. (2007).
\emph{Lehrbuch Theorie und Methoden der Skalierung}.
Bern: Hogrefe.
}


\examples{
\dontrun{
#############################################################################
# EXAMPLE 07a: Dataset Gefechtsangst
#############################################################################

data(data.bs07a)
dat <- data.bs07a
items <- grep( "GF" , colnames(dat)  , value=TRUE )

#************************
# Model 1: Rasch model
mod1 <- TAM::tam.mml(dat[,items] )
summary(mod1)
IRT.WrightMap(mod1)

#************************
# Model 2: 2PL model
mod2 <- TAM::tam.mml.2pl(dat[,items] )
summary(mod2)

#************************
# Model 3: Latent class analysis (LCA) with two classes
tammodel <- "
ANALYSIS:
  TYPE=LCA;
  NCLASSES(2)
  NSTARTS(5,10)
LAVAAN MODEL:
  F =~ GF1__GF9
  "
mod3 <- TAM::tamaan( tammodel , dat )
summary(mod3)

#************************
# Model 4: LCA with three classes
tammodel <- "
ANALYSIS:
  TYPE=LCA;
  NCLASSES(3)
  NSTARTS(5,10)
LAVAAN MODEL:
  F =~ GF1__GF9
  "
mod4 <- TAM::tamaan( tammodel , dat )
summary(mod4)

#************************
# Model 5: Located latent class model (LOCLCA) with two classes
tammodel <- "
ANALYSIS:
  TYPE=LOCLCA;
  NCLASSES(2)
  NSTARTS(5,10)
LAVAAN MODEL:
  F =~ GF1__GF9
  "
mod5 <- TAM::tamaan( tammodel , dat )
summary(mod5)

#************************
# Model 6: Located latent class model with three classes
tammodel <- "
ANALYSIS:
  TYPE=LOCLCA;
  NCLASSES(3)
  NSTARTS(5,10)
LAVAAN MODEL:
  F =~ GF1__GF9
  "
mod6 <- TAM::tamaan( tammodel , dat )
summary(mod6)

#************************
# Model 7: Probabilistic Guttman model
mod7 <- sirt::prob.guttman( dat[,items] )
summary(mod7)

#-- model comparison
IRT.compareModels( mod1, mod2 , mod3 , mod4 , mod5 , mod6 , mod7 )
}
}
\keyword{datasets}
