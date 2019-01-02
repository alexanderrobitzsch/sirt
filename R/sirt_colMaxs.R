## File Name: sirt_colMaxs.R
## File Version: 0.06

sirt_colMaxs <- function(x, na.rm=TRUE)
{
    res <- apply(x, 2, max, na.rm=na.rm)
    return(res)
}

sirt_colMax <- sirt_colMaxs
