## File Name: automatic.recode.R
## File Version: 1.194

#*** automatic recoding of a dataset
automatic.recode <- function( data, exclude=NULL, pstart.min=.6,
    allocate=200, maxiter=20, progress=TRUE)
{
    # data table with frequencies and discriminations
    item.stat <- data.frame( item=colnames(data), old.key=-99, new.key=NA,
                                p=NA, discrim=NA )
    items <- colnames(data)
    # compute frequencies
    fstart <- TAM::tam.ctt3( data, allocate=allocate, progress=FALSE)
    I <- ncol(data)
    prbar <- floor( 10 * ( 1L:I ) / (I+1) )
    prbar <- c(1,which( diff(prbar)==1 )  )

    fstart1 <- fstart
    fstart1 <- fstart1[ ! ( fstart$Categ %in% exclude ), ]
    fstart1 <- fstart1[ fstart1$RelFreq > pstart.min, ]
    fstart1 <- fstart1[ order( fstart1$itemno + fstart1$RelFreq, decreasing=TRUE), ]
    fstart1 <- fstart1[ ! duplicated( fstart1$itemno), ]
    item.stat[ fstart1$itemno, 'new.key' ] <- fstart1[, 'Categ' ]
    item.stat[ fstart1$itemno, 'p' ] <- fstart1[, 'RelFreq' ]
    iter <- 1
    changed.keys <- I

    #--- start iterations
    while ( ( iter < maxiter ) & ( changed.keys > 0 ) ){
        # select items which should be recoded
        is1 <- paste( item.stat[ ! is.na( item.stat$new.key ), 'item' ] )

        # scoring responses
        data.raw <- data[, is1 ]
        item.stat$key <- item.stat$new.key
        # scoring data
        data.scored <- data.recode.sirt( data.raw, keys=item.stat )

        score <- rowMeans( data.scored, na.rm=TRUE )
        # compute frequencies in iterations
    #    if (progress){ cat( '   apply tam.ctt3\n') }
        fiter <- TAM::tam.ctt3( data, allocate=allocate, wlescore=score, progress=FALSE)
    #    if (progress){ cat( '   |') }
        item.stat$old.key <- item.stat$new.key

        fstart1 <- fiter
        fstart1 <- fstart1[ ! ( fstart$Categ %in% exclude ), ]
        fstart1 <- fstart1[ order( fstart1$itemno + fstart1$rpb.WLE, decreasing=TRUE), ]
        fstart1 <- fstart1[ ! duplicated( fstart1$itemno), ]

        item.stat[ fstart1$itemno, 'new.key' ] <- fstart1[, 'Categ' ]
        item.stat[ fstart1$itemno, 'p' ] <- fstart1[, 'RelFreq' ]
        item.stat[ fstart1$itemno, 'discrim' ] <- fstart1[, 'rpb.WLE' ]

        changed.keys <- sum( paste(item.stat$new.key) !=paste(item.stat$key) )
        if (progress){
            cat( paste0( 'Iteration ', iter,     ' | ',
                            'Changed ', changed.keys, ' Keys\n' ) )
            utils::flush.console()
        }
        iter <- iter+1
    }
    #********* end iterations
    item.stat$old.key <- item.stat$new.key <- NULL
    res <- list( item.stat=item.stat, data.scored=data.scored, categ.stats=fiter)
    return(res)
}
