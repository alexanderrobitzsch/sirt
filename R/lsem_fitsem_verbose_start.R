## File Name: lsem_fitsem_verbose_start.R
## File Version: 0.05
## File Last Change: 2023-03-11

lsem_fitsem_verbose_start <- function(G, verbose)
{
    pr <- NULL
    if (verbose){
        cat( '** Fit lavaan model\n')
        G1 <- min(G,10)
        pr <- round( seq(1, G, len=G1) )
        cat('|')
        cat( paste0( rep('*',G1), collapse='') )
        cat('|\n')
        cat('|')
    }
    return(pr)
}
