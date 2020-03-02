## File Name: xxirt_irf_lca_init_par.R
## File Version: 0.03


xxirt_irf_lca_init_par <- function(K, ncat, random_sd=0, rg_val=1.5)
{
    v1 <- NULL
    for (cc in 2L:ncat){
        par <- rg_val*seq( -1, 1, length=K )
        names(par) <- paste0("b",1:K, "_Cat",cc-1)
        par <- par + stats::rnorm(K, sd=random_sd)
        v1 <- c(v1, par)
    }
    if (ncat==2){
        names(v1) <- paste0("b",1:K)
    }
    return(v1)
}
