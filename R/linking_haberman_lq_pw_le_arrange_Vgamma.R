## File Name: linking_haberman_lq_pw_le_arrange_Vgamma.R
## File Version: 0.04


linking_haberman_lq_pw_le_arrange_Vgamma <- function(vcov_list, par_gamma, I, G,
        ind_items, ind_studies)
{

    NPG <- 2*I*G
    Vgamma <- matrix(0, nrow=NPG, ncol=NPG,
                    dimnames=list(names(par_gamma), names(par_gamma) ) )
    if (!is.null(vcov_list)){
        for (gg in 1L:G){
            items_gg <- ind_items[ ind_studies==gg ]
            Igg <- length(items_gg)
            items_gg <- rep(items_gg, each=2)
            ind2 <- 2*I*(gg-1) + 2*(items_gg-1)+rep(1L:2, Igg)
            Vgamma[ind2, ind2] <- vcov_list[[gg]]
        }
    }
    return(Vgamma)
}
