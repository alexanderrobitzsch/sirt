## File Name: lsem_fitsem_verbose_progress.R
## File Version: 0.03

lsem_fitsem_verbose_progress <- function(gg, G, pr, verbose)
{
    if (verbose){
        if (gg %in% pr ){
            cat('-')
        }
        if (gg>G){
            cat('|\n')
        }
        utils::flush.console()
    }
}
