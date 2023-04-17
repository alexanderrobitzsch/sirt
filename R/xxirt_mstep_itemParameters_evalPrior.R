## File Name: xxirt_mstep_itemParameters_evalPrior.R
## File Version: 0.153


#*** evaluate prior in M-step
xxirt_mstep_itemParameters_evalPrior <- function(partable, h=0)
{
    eps <- 1E-300
    partable1 <- partable[ partable$parfree==1, ]
    NP <- nrow(partable1)
    NP1 <- max( NP, 1)
    #*** evaluate prior distributions in partable
    pen <- rep(0,NP1)
    if (NP>0){
        for (pp in 1:NP){
            if ( ! is.na( partable1[pp,'prior'] ) ) {
                prior_pp <- partable1[pp,'prior']
                val <- partable1[pp,'value'] + h
                prior_args_pp <- list( val, partable1[pp,'prior_par1'],
                                                partable1[pp,'prior_par2'] )

                prior_val <- do.call( prior_pp, prior_args_pp )
                pen[pp] <- - log( prior_val + eps )
            }  # end if prior
        }    # end pp

    }  # end if NP > 0
    return(pen)
}
