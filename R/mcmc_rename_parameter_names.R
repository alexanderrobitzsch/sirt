## File Name: mcmc_rename_parameter_names.R
## File Version: 0.121

mcmc_rename_parameter_names <- function( vec, orig, trans)
{
    NO <- length(orig)
    for (oo in 1L:NO){
        trans_oo <- mcmc_rename_helper( string=trans[oo] )
        trans_oo <- gsub( ' ', '', trans_oo)
        vec <- gsub( orig[oo], trans_oo, vec, fixed=TRUE)
    }
    return(vec)
}
