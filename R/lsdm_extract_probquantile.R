## File Name: lsdm_extract_probquantile.R
## File Version: 0.14


# auxiliary function for extracting quantiles of curves
lsdm_extract_probquantile <- function( vec, theta, quant )
{
    x2 <- theta[ vec>=quant ][1]
    x1 <- sort( theta[ vec<quant ], decreasing=TRUE)[1]
    value <- - Inf
    if ( (1-is.na(x2)) * (1-is.na(x1))==1 ){
        y1 <- vec[theta==x1]
        y2 <- vec[theta==x2]
        value <- x1 + ( quant - y1 ) * ( x2 - x1 ) / ( y2 - y1 )
    }
    if ( is.na(x2) ){
        value <- Inf
    }
    return(value)
}
