## File Name: dmlavaan_adjust_numdiff_h.R
## File Version: 0.03

dmlavaan_adjust_numdiff_h <- function(h, val)
{
    # h1 <- ifelse(abs(val)>1, h*abs(val), h )
    h1 <- h
    return(h1)
}
