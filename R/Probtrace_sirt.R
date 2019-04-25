## File Name: Probtrace_sirt.R
## File Version: 0.05


# auxiliary function for mirt package
# trace function for all items
Probtrace_sirt <- function(items, Theta)
{
    TAM::require_namespace_msg("mirt")
    traces <- lapply(items, mirt::probtrace, Theta=Theta)
    ret <- do.call(cbind, traces)
    return(ret)
}
