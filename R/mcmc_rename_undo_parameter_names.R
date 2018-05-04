## File Name: mcmc_rename_undo_parameter_names.R
## File Version: 0.01

mcmc_rename_undo_parameter_names <- function( vec, orig, trans)
{
    NO <- length(orig)
    for (oo in 1:NO){
        trans_oo <- mcmc_rename_helper( string=trans[oo] )        
        vec <- gsub( trans_oo, orig[oo], vec , fixed=TRUE)
    }
    return(vec)
}
