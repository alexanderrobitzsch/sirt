## File Name: regpolca_postproc_count_regularized_parameters.R
## File Version: 0.28
## File Last Change: 2020-07-04

regpolca_postproc_count_regularized_parameters <- function(item, set_equal,
        lca_dich, probs_items, nclasses, ncats)
{
    item1_index <- NULL
    n_reg <- 0
    if (lca_dich){
        ind <- c(1)
        item1 <- item[,-ind]
        TP <- ncol(item1)
        I <- nrow(item)
        item1_index <- 0*item1

        for (ii in 1:I){
            index_ii <- item1_index[ii,]
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
    if ((! lca_dich) & (set_equal>0)){
        items <- unique(paste(item$item))
        I <- length(items)
        n_reg <- 0
        lem1 <- c(1,2)
        for (ii in 1:I){
            item_ii <- items[ii]
            ind_ii <- which( item[,"item"]==item_ii )
            item_tab_ii <- item[ind_ii,-lem1]
            item_tab_ii[-c(ind_ii[1]),] <- set_equal*round( item_tab_ii[-c(ind_ii[1]),] / set_equal, 0 )
            item_tab_ii[1,] <- 1 - colSums(item_tab_ii[-1,])
            item[ind_ii, -c(lem1) ] <- item_tab_ii
            LU <- length( unique( unlist(item_tab_ii[-1,] )))
            H <- (nrow(item_tab_ii)-1)*ncol(item_tab_ii)
            n_reg <- n_reg + H - LU
            for (hh in 1:ncats[ii]){
                v1 <- as.numeric(item_tab_ii[hh,])
                probs_items[ii,hh,] <- v1
            }
        }
    }

    res <- list(item1_index=item1_index, n_reg=n_reg, item=item,
                probs_items=probs_items)
    return(res)
}
