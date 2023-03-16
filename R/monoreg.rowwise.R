## File Name: monoreg.rowwise.R
## File Version: 0.11
## File Last Change: 2018-12-30


# monotone regression for all rows in a matrix
monoreg.rowwise <- function(yM, wM)
{
    yM <- as.matrix(yM)
    wM <- as.matrix(wM)
    res <- sirt_rcpp_monoreg_rowwise( YM=yM, WM=wM )
    return(res)
}

