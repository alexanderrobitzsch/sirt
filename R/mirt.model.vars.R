## File Name: mirt.model.vars.R
## File Version: 0.259


#-- Input is a mirt syntax specified with variable names instead of numbers
mirt.model.vars <- function( mirtsyn2, data=NULL, items=NULL )
{
    TAM::require_namespace_msg('mirt')
    # mirt syntax splitted
    mirtsyn2 <- gsub( ';', '\n', mirtsyn2 )
    msyn0 <- strsplit( mirtsyn2, c(' ') )[[1]]
    syn <- msyn0
    if ( !is.null(data) ){
        items <- colnames(data)
    }
    if ( is.null(items) ){
        stop('Provide dataset or items.\n')
    }
    # admissible strings
    vecstr <- c( '\n', '\\(', '=', '-', ',', '\\)' )

    #--- process syntax
    for (vv in vecstr){
        syn <- split_syn_string( syn, vv )
    }
    # postprocess syntax
    syn <- syn[ syn !='' ]
    syn[ syn=='\\(' ] <- '('
    syn[ syn=='\\)' ] <- ')'
    # substitute variables by numbers
    VV <- length(items)
    useditems <- NULL
    for (vv in 1L:VV){
        ind <- which( syn==items[vv] )
        if (length(ind) > 0 ){
            syn[ ind ] <- vv
            useditems <- c( useditems, items[vv] )
        }
    }
    syn <- paste0( syn, collapse='')
    mirtmodel <- mirt::mirt.model(syn)
    mirtmodel <- list(model=mirtmodel, syntax=syn )
    return(mirtmodel)
}


