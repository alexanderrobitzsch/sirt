## File Name: rm_sdt_fill_init_partable.R
## File Version: 0.09


rm_sdt_fill_init_partable <- function(partable, par_index, tau.item, a.item,
    c.rater, d.rater, type)
{
    #--- collect parameters
    if( type=='item' ){
        pars <- list( tau.item=tau.item, a.item=a.item)
    } else {
        pars <- list( c.rater=c.rater, d.rater=d.rater)
    }
    #--- inits
    for (pp in names(pars)){
        par_pp <- par_index[[pp]]
        if ( length(par_pp) > 0 ){
            partable[ par_pp, "value" ] <- as.vector( pars[[ pp ]] )
        }
    }
    #-- output
    return(partable)
}
