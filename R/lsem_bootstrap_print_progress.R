## File Name: lsem_bootstrap_print_progress.R
## File Version: 0.08

lsem_bootstrap_print_progress <- function(rr, verbose, R)
{
    if (verbose){
        cat(rr, " ")
        if (rr %% 10==0){ cat("\n")}
        if (rr==R){ cat("\n") }
        utils::flush.console()
    }
}
