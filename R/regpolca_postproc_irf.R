## File Name: regpolca_postproc_irf.R
## File Version: 0.04


regpolca_postproc_irf <- function(probs_items, dat, lca_dich)
{
    K <- dim(probs_items)[3]
    NC <- dim(probs_items)[2]
    I <- dim(probs_items)[1]
    ci <- paste0("Cl",1:K)
    if (lca_dich){
        item <- probs_items[,2,]
        colnames(item) <- ci
        item <- data.frame( item=colnames(dat), item)
    } else {
        item <- NULL
        for (cc in 1:NC){
            dfr1 <- probs_items[,cc,]
            colnames(dfr1) <- ci
            dfr1 <- data.frame(item=1:I, cat=paste0("Cat",cc-1), dfr1)
            item <- rbind(item, dfr1)
        }
        item <- item[ order(item$item), ]
        item$item <- colnames(dat)[ item$item ]
    }
    return(item)
}
