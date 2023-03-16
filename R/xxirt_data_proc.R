## File Name: xxirt_data_proc.R
## File Version: 0.213
## File Last Change: 2023-02-15

#-- data processing xxirt
xxirt_data_proc <- function(dat, group=NULL, weights=NULL )
{
    ncat <- apply( dat, 2, max, na.rm=TRUE ) + 1
    I <- ncol(dat)
    items <- colnames(dat)
    N <- nrow(dat)
    if ( is.null(group) ){
        group <- rep(1,N)
    }
    if ( is.null(weights) ){
        weights <- rep(1,N)
    }
    W <- sum(weights)
    group0 <- group
    groups_unique <- sort( unique( group ) )
    G <- length(groups_unique)
    group <- match( group0, groups_unique )
    maxK <- max(ncat)
    #*** group_index
    group_index <- as.list( 1:G )
    for (gg in 1:G){
        group_index[[gg]] <- which( group==gg )
    }
    #*** data with response indices
    dat_na <- is.na(dat)
    dat_resp_bool <- ! dat_na
    dat_resp <- 1 - dat_na
    resp_index <- as.list( 1:I )
    for ( ii in 1:I){
        resp_index[[ii]] <- which( dat_resp[,ii]==1 )
    }
    dat1 <- as.matrix(dat)
    dat1[ dat_na ] <- 0
    #*** output
    res <- list( N=N, I=I, group=group, items=items,
                    group0=group0, G=G, groups_unique=groups_unique,
                    maxK=maxK, ncat=ncat, weights=weights,
                    group_index=group_index, dat_resp=dat_resp,
                    resp_index=resp_index, dat1=dat1, dat_resp_bool=dat_resp_bool,
                    W=W)
    return(res)
}

