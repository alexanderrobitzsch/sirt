## File Name: gom_em_prepare_lambda_index.R
## File Version: 0.12

gom_em_prepare_lambda_index <- function(lambda.index, I, K, items)
{
    if (is.null(lambda.index)){
        lambda.index <- matrix(1:(I*K), nrow=I, ncol=K)
    }
    li <- as.vector(lambda.index)
    lambda.index <- matrix( match(li, unique(as.vector(li))), nrow=I, ncol=K)
    lambda_partable <- data.frame(par_index=as.vector(lambda.index),
                        item=rep(1:I, K), class=rep(1:K, each=I) )
    lambda_partable$free <- 1
    lambda_partable[ duplicated(lambda_partable$par_index), "free"] <- 0
    lambda_partable$block <- 0

    #- check whether there are constraints among item parameters
    item_con <- any( unlist( apply( lambda.index, 2, FUN=function(ll){ length(unique(ll))<I } )))
    lambda_partable$item_con <- 1*item_con

    lambda_index_blocks    <- list()
    bb <- 0
    for (kk in 1:K){
        pkk <- lambda_partable[ (lambda_partable$class==kk ) & (lambda_partable$free==1), "par_index"]
        if (length(pkk)>0){
            bb <- bb+1
            lambda_index_blocks[[bb]] <- pkk
            lambda_partable[ lambda_partable$par_index %in% pkk, "block" ] <- bb
        }
    }
    lambda_partable$block <- lambda_partable$block * (lambda_partable$free > 0)
    if (item_con){
        lambda_partable$block <- lambda_partable$par_index * (lambda_partable$free > 0)
        lambda_index_blocks <- NULL
    }

    #-- parameter names
    lambda_partable$par_name <- paste0( rep(items, K), "_Cl", rep(1:K, each=I))

    #-- output
    res <- list(lambda.index=lambda.index, lambda_partable=lambda_partable,
                lambda_index_blocks=lambda_index_blocks)
    return(res)
}
