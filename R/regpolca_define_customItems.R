## File Name: regpolca_define_customItems.R
## File Version: 0.06

regpolca_define_customItems <- function(ncats, K, dat, par_item_max)
{
    cats_u <- unique(ncats)
    customItems <- list(1:length(cats_u))
    for (ncat in cats_u){
        par <- xxirt_irf_lca_init_par(K=K, ncat=ncat, random_sd=0)
        item_lc <- xxirt_createDiscItem( name=paste0('LC',ncat), par=par,
                     est=rep(TRUE,K*(ncat-1)), P=xxirt_irf_lca )
        ii <- which(ncat==cats_u)
        customItems[[ii]] <- item_lc
    }
    itemtype <- paste0( 'LC', ncats )
    partable <- xxirt_createParTable( dat=dat, itemtype=itemtype,
                    customItems=customItems)
    partable$lower <- -par_item_max
    partable$upper <- par_item_max

    #- output
    res <- list(customItems=customItems, partable=partable,
                    itemtype=itemtype)
    return(res)
}

