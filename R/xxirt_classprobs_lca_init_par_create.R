## File Name: xxirt_classprobs_lca_init_par_create.R
## File Version: 0.04
## File Last Change: 2020-02-24


xxirt_classprobs_lca_init_par_create <- function(K, random_sd=0)
{
    par_Theta <- rep(0, K-1)
    par_Theta <- par_Theta + stats::rnorm(K-1, sd=random_sd)
    return(par_Theta)
}
