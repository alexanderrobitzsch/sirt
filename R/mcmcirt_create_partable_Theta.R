## File Name: mcmcirt_create_partable_Theta.R
## File Version: 0.13
## File Last Change: 2022-10-23


mcmcirt_create_partable_Theta <- function(par, est, prior=NULL,
        prior_par1=NULL, prior_par2=NULL, sd_proposal=NULL)
{
    NP <- length(par)
    par_names <- names(par)
    dfr <- NULL
    for (pp in 1L:NP){
        name_pp <- par_names[pp]
        dfr1 <- data.frame(index=pp, par=name_pp, start=par[pp])
        dfr1$est <- est[name_pp]
        #** prior distribution
        dfr1$prior <- prior[name_pp]
        dfr1$prior_par1 <- prior_par1[name_pp]
        dfr1$prior_par2 <- prior_par2[name_pp]
        dfr1$sd_proposal <- sd_proposal[name_pp]
        dfr <- rbind(dfr, dfr1)
    }
    dfr$parindex <- cumsum(est)
    dfr$value <- dfr$start
    dfr$sampled <- 0
    dfr$accepted <- 0
    dfr$EAP <- 0
    return(dfr)
}
