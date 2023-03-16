## File Name: monoreg.colwise.R
## File Version: 0.13


# monotone regression for all columns in a matrix
monoreg.colwise <- function(yM, wM)
{
    yM <- as.matrix(t(yM))
    wM <- as.matrix(t(wM))
    res <- sirt_rcpp_monoreg_rowwise( YM=yM, WM=wM )
    return(t(res))
}
