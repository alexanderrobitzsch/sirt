## File Name: lsem_aggregate_statistics.R
## File Version: 0.02


lsem_aggregate_statistics <- function(x)
{
    Nimp <- length(x)
    w <- 0
    for (ii in 1L:Nimp){
        w <- w + x[[ii]]
    }
    w <- w / Nimp
    return(w)
}
