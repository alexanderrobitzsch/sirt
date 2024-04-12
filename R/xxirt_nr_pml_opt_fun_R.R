## File Name: xxirt_nr_pml_opt_fun_R.R
## File Version: 0.06


xxirt_nr_pml_opt_fun_R <- function(prior_Theta, probs_items, freq1, freq2, W1,
                                W2_long, G, K, I, TP, NI2, eps )
{
    val <- 0

    #** first-order frequencies
    for (gg in 1L:G){
        for (ii in 1L:I){
            t1 <- 0
            for (hh in 1L:K){
                p1 <- sum( prior_Theta[,gg] * probs_items[ii,hh,] )
                t1 <- t1 + freq1[ii,hh,gg]*log(p1+eps)
            }
            val <- val + W1[ii]*t1
        }
    }

    #*** second-order frequencies
    for (gg in 1L:G){
        for (nn in 1L:NI2){
            ii1 <- W2_long[nn,'item1']
            ii2 <- W2_long[nn,'item2']
            w2ii <- W2_long[nn,'w2']
            t1 <- 0
            for (hh in 1L:K){
                for (kk in 1L:K){
                    p1 <- sum( prior_Theta[,gg] * probs_items[ii1,hh,] *
                                                probs_items[ii2,kk,])
                    t1 <- t1 + freq2[nn,hh,kk,gg]*log(p1+eps)
                }
            }
            val <- val + w2ii*t1
        }
    }

    #- output
    res <- list(vcal=val)
    return(res)
}
