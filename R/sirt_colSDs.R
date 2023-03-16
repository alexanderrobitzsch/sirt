## File Name: sirt_colSDs.R
## File Version: 0.04
## File Last Change: 2018-12-30

sirt_colSDs <- function(x, na.rm=TRUE)
{
    res <- apply(x, 2, stats::sd, na.rm=na.rm)
    return(res)
}

sirt_colSD <- sirt_colSDs
