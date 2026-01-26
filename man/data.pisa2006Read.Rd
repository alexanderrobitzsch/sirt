%% File Name: data.pisa2006Read.Rd
%% File Version: 0.032

\name{data.pisa2006Read}
\alias{data.pisa2006Read}
\docType{data}
\title{
Dataset PISA 2006 Reading
}
\description{
This is the PISA 2006 dataset of item responses from 26 selected OECD countries.
Item responses from the reading domain are selected.
}
\usage{
data(data.pisa2006Read)
}


\format{
The dataset includes item responses, student weights and their respective
replicate weights.

The format is

\code{'data.frame':   110236 obs. of  116 variables:} \cr
\code{ $ index   : int  1 2 3 4 5 6 7 8 9 10 ...} \cr
\code{ $ cnt     : Factor w/ 26 levels "AUS","AUT","BEL",..: 1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{ $ IDSTUD  : int  603600001 603600002 603600006 603600008 603600009 603600010 603600011 603600012 603600013 603600014 ...} \cr
\code{ $ COUNTRY : int  36 36 36 36 36 36 36 36 36 36 ...} \cr
\code{ $ SCHOOLID: int  1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{ $ STIDSTD : int  1 2 6 8 9 10 11 12 13 14 ...} \cr
\code{ $ BOOKID  : int  2 13 11 2 6 13 12 7 12 11 ...} \cr
\code{ $ wgt     : num  0.239 0.239 0.218 0.218 0.218 ...} \cr
\code{ $ R055Q01 : int  NA 1 1 NA 1 1 NA NA NA 1 ...} \cr
\code{ $ R055Q02 : int  NA 0 0 NA 1 1 NA NA NA 0 ...} \cr
\code{ $ R055Q03 : int  NA 1 1 NA 1 1 NA NA NA 1 ...} \cr
\code{ $ R055Q05 : int  NA 0 1 NA 1 1 NA NA NA 1 ...} \cr
\code{[...]} \cr
\code{ $ R227Q03 : int  NA 0 0 NA 1 1 NA NA NA 0 ...} \cr
\code{ $ R227Q06 : int  NA 1 1 NA 1 1 NA NA NA 0 ...} \cr
\code{ $ fac_rw1 : num  1.47 0.49 1.49 1.49 1.49 1.49 0.49 0.5 0.5 0.49 ...} \cr
\code{ $ fac_rw2 : num  1.67 0.56 1.57 1.57 0.52 0.52 0.56 0.52 1.57 0.56 ...} \cr
\code{ $ fac_rw3 : num  1.46 1.46 0.5 1.49 0.5 1.49 1.46 1.49 1.49 1.46 ...} \cr
\code{[...]} \cr
}



%% \keyword{datasets}
