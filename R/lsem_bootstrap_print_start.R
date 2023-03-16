## File Name: lsem_bootstrap_print_start.R
## File Version: 0.02
## File Last Change: 2023-03-15

lsem_bootstrap_print_start <- function(verbose)
{
    if (verbose){
        cat('Compute LSEM for bootstrap samples:\n')
        utils::flush.console()
    }
}
