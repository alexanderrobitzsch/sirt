## File Name: regpolca_penalty_fun.R
## File Version: 0.05


regpolca_penalty_fun <- function(x, regular_grouped, I, partable,
    combis_classes_list, regular_lam, eps, penalty_used, Ni,
    combis_categ_list, fuse_categories, K)
{
    pen <- 0
    if (regular_grouped=='none'){
        #* fused probabilities among classes
        for (ii in 1L:I){
            x_ii <- x[ partable$itemnr==ii ]
            a1 <- regpolca_penalty_fun_value_nongrouped(x_ii=x_ii,
                                combis_ii=combis_classes_list[[ii]],
                                regular_lam=regular_lam[1],
                                eps=eps, penalty_used=penalty_used)
            pen <- pen + Ni[ii]*sum(a1)
        }
        #* fixed probabilities among categories
        for (ii in 1L:I){
            if (fuse_categories[ii]){
                x_ii <- x[ partable$itemnr==ii ]
                a1 <- regpolca_penalty_fun_value_nongrouped(x_ii=x_ii,
                            combis_ii=combis_categ_list[[ii]],
                            regular_lam=regular_lam[2],
                            eps=eps, penalty_used=penalty_used)
                pen <- pen + Ni[ii]*sum(a1)
            }
        }
    } else {
        for (ii in 1L:I){
            x_ii <- x[ partable$itemnr==ii ]
            combis_ii <- combis_categ_list[[ii]]
            nc <- nrow(combis_ii)
            a1 <- 0
            if (regular_grouped=='class'){
                K1 <- max(combis_ii$cat_pair)
                for (kk in 1:K1){
                    cii_kk <- combis_ii[ seq(kk,nc,by=nc/K), ]
                    a2 <- regpolca_penalty_fun_value_grouped(x_ii=x_ii,
                                combis_ii=cii_kk, regular_lam=regular_lam[1],
                                eps=eps, penalty_used=penalty_used)
                    a1 <- a1 + a2
                }
            }
            if (regular_grouped=='categ'){
                for (kk in 1:K){
                    cii_kk <- combis_ii[ combis_ii$class==kk, ]
                    a2 <- regpolca_penalty_fun_value_grouped(x_ii=x_ii,
                                combis_ii=cii_kk, regular_lam=regular_lam[1],
                                eps=eps, penalty_used=penalty_used)
                    a1 <- a1 + a2
                }
            }
            pen <- pen + Ni[ii]*sum(a1)
        }
    }
    return(pen)
}
