## File Name: mgsem_vech.R
## File Version: 0.01

mgsem_vech <- function(x)
{
    res <- x[ ! upper.tri(x)]
    return(res)
}
