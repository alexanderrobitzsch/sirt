## File Name: xxirt_createParTable.R
## File Version: 0.254
## File Last Change: 2023-03-08

#*** create parameter table
xxirt_createParTable <- function( dat, itemtype, customItems=NULL )
{
    I <- ncol(dat)
    ncat1 <- apply( dat, 2, max, na.rm=TRUE ) + 1
    items <- colnames(dat)
    if ( length(itemtype)==1 ){
        itemtype <- rep( itemtype, I )
    }
    dfr <- NULL
    CI <- length(customItems)
    for (ii in 1:I){
        type_ii <- itemtype[ii]
        item_ii <- NULL
        for ( vv in 1:CI){
            ci_ii <- customItems[[vv]]
            if ( ci_ii$name==type_ii ){
                item_ii <- ci_ii
            }
        }
        if ( is.null(item_ii) ){
            stop( paste0( 'Item type ', type_ii, ' not found!') )
        }
        NP <- length( item_ii$par )
        dfr1 <- data.frame( 'item'=rep( items[ii], NP ) )
        dfr1$itemnr <- ii
        dfr1$ncat <- ncat1[ii]
        dfr1$class <- class(item_ii)
        dfr1$type <- item_ii$name
        dfr1$parname <- names(item_ii$par)
        dfr1$value <- item_ii$par
        dfr1$est <- item_ii$est
        dfr1$lower <- item_ii$lower
        dfr1$upper <- item_ii$upper
        dfr1$prior <- NA
        dfr1$prior_par1 <- NA
        dfr1$prior_par2 <- NA
        if ( ! is.null( item_ii$prior ) ){
            item_ii_prior <- names(item_ii$prior)
            ind_ii <- match( item_ii_prior, dfr1$parname )
            dfr1[ ind_ii, 'prior' ] <- item_ii$prior
            dfr1[ ind_ii, 'prior_par1' ] <- item_ii$prior_par1
            dfr1[ ind_ii, 'prior_par2' ] <- item_ii$prior_par2
        }
        dfr <- rbind( dfr, dfr1 )
    }
    #**** create parameter indices
    NP <- nrow(dfr)
    dfr$rowindex <- 1:NP
    # parameter index
    dfr$parindex <- cumsum( dfr$est )
    #*** parameter label
    dfr$parlabel <- paste0( dfr$item, '_', dfr$parname )
    attr(dfr, 'ncat' ) <- ncat1
    attr(dfr, 'items' ) <- items
    return(dfr)
}
