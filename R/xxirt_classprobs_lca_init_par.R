## File Name: xxirt_classprobs_lca_init_par.R
## File Version: 0.091
## File Last Change: 2023-03-08


xxirt_classprobs_lca_init_par <- function(K, G, random_sd=0, par_Theta_init=NULL)
{
    if (G==1){
        par_Theta <- xxirt_classprobs_lca_init_par_create(K=K, random_sd=random_sd)
        names(par_Theta) <- paste0('pi',1:(K-1) )
        if ( !is.null(par_Theta_init) ){
            par_Theta <- par_Theta_init
        }
    } else {
        v1 <- NULL
        for (gg in 1L:G){
            par_Theta <- xxirt_classprobs_lca_init_par_create(K=K, random_sd=random_sd)
            names(par_Theta) <- paste0('pi',1:(K-1), '_Gr', gg)
            v1 <- c(v1, par_Theta)
        }
        par_Theta <- v1
    }
    return(par_Theta)
}
