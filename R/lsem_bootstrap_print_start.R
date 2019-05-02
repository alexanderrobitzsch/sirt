## File Name: lsem_bootstrap_print_start.R
## File Version: 0.01

lsem_bootstrap_print_start <- function(verbose)
{
    if (verbose){
        cat("Compute LSEM for bootstrap samples:\n")
        utils::flush.console()
    }
}
