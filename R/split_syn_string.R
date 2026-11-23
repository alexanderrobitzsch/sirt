## File Name: split_syn_string.R
## File Version: 0.01


#-- cleans syntax in a vector from strings vv
split_syn_string <- function( syn, vv )
{
    syn <- as.list(syn )
    syn.vv <- grep( vv, syn )
    LL <- length(syn.vv)
    if (LL>0){
        for (ii in 1:LL){
            ll <- syn.vv[ii]
            syn.ll <- syn[[ll]]
            syn[[ll]] <- split_conc( syn.ll, vv )
        }
    }
    syn <- unlist(syn)
    return(syn)
}
