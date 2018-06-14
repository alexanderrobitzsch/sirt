## File Name: mcmc_derivedPars.R
## File Version: 0.23

#####################################################
# derived parameters for objects of class mcmc
mcmc_derivedPars <- function( mcmcobj, derivedPars )
{
    mcmcobj <- mcmc_extract_samples_first_chain(mcmcobj=mcmcobj)
    NP <- length(derivedPars)
    data <- as.data.frame( mcmcobj )
    #-- renaming
    res <- mcmc_rename_define_symbols()
    orig <- res$orig
    trans <- res$trans
    colnames(data) <- mcmc_rename_parameter_names( vec=colnames(data), orig=orig, trans=trans)
    for (pp in 1:NP){
        der_pp <- mcmc_rename_parameter_names( vec=paste0(derivedPars[[pp]]), orig=orig, trans=trans)
        form_pp <- mcmc_as_formula(der_pp)
        data_pp <- stats::model.matrix( form_pp , data )
        if (ncol(data_pp) > 1){
            data_pp <- data_pp[,-1]
        }
        data <- as.data.frame( cbind( data, data_pp) )
        colnames(data)[ ncol(data)] <- names(derivedPars)[pp]
    }
    colnames(data) <- mcmc_rename_undo_parameter_names( vec=colnames(data), orig=orig, trans=trans)
    a1 <- attr( mcmcobj , "mcpar")
    res <- coda::mcmc(     data= data ,  start = a1[1] , end = a1[2],
                    thin = a1[3] )
    return(res)
}
