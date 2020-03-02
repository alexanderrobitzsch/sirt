## File Name: regpolca_postproc_count_regularized_parameters.R
## File Version: 0.08

regpolca_postproc_count_regularized_parameters <- function(item, set_equal,
        lca_dich, probs_items)
{

    if (lca_dich){
        ind <- c(1)
        item1 <- item[,-ind]
        TP <- ncol(item1)
        I <- nrow(item)
        item1_index <- 0*item1
        index_ii <- item1_index[ii,]
        for (ii in 1:I){
            item_ii <- item1[ii,]
            NC <- length(item_ii)
            index_ii <- 1:NC
            for (aa in 1:(NC-1)){
                val_aa <- item_ii[aa]
                for (bb in (aa+1):NC){
                    val_bb <- item_ii[bb]
                    if (abs(val_aa-val_bb)<set_equal){
                        index_ii[bb] <- index_ii[aa]
                    }
                }
            }
            ind_ii <- match( index_ii, unique(index_ii) )
            item1_index[ii,] <- ind_ii
            w1 <- t(item_ii)
            a1 <- stats::aggregate(w1[,1], list(ind_ii), mean)
            probs_mod <- a1[ind_ii,2]
            item[ii, -ind] <- probs_mod
            probs_items[ii,2,] <- probs_mod
            probs_items[ii,1,] <- 1 - probs_mod
        }  # end ii
        H <- I*TP
        n_reg <- H - sum( apply(item1_index, 1, max))
    }
    res <- list(item1_index=item1_index, n_reg=n_reg, item=item,
                probs_items=probs_items)
    return(res)
}
