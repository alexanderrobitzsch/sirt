## File Name: ccov.np.R
## File Version: 1.221


#---- nonparametric estimation of conditional covariance
ccov.np <- function( data, score, bwscale=1.1, thetagrid=seq( -3,3,len=200),
        progress=TRUE, scale_score=TRUE, adjust_thetagrid=TRUE, smooth=TRUE,
        use_sum_score=FALSE, bias_corr=TRUE)
{
    # number of Items I
    I <- ncol(data)
    if (use_sum_score){
        smooth <- FALSE
        score <- rowSums(data)
    }
    # z-standardization of score
    if (scale_score){
        score <- scale(score)[,1]
    } else {
        bwscale <- bwscale*stats::sd(score, na.rm=TRUE)
    }
    # adjust thetagrid
    if (adjust_thetagrid & smooth){
        score_rg <- range(score)
        thetagrid <- seq(score_rg[1], score_rg[2], length=length(thetagrid))
    }
    if (!smooth){
        score <- round(score,3)
        thetagrid <- sort(unique(score))
    }
    # matrix of item response functions
    if (progress){
        cat('Pairwise Estimation of Conditional Covariances\n' )
        cat('...........................................................\n' )
        cat('Nonparametric ICC estimation \n ' )
    }
    icc_items <- matrix( 0, length(thetagrid), I )
    if ( I >=20 ){
        display <- seq( 1, I, floor( I/20 ) )[ 2:20 ]
    } else {
        display <- 20
    }
    i <- 1
    for ( ii in 1:I ){
        obs_ii <- ! is.na( data[,ii] )
        x <- score[ obs_ii ]
        y <- data[ obs_ii, ii ]
        icc_items[,ii] <- ccov_np_regression(x=x, y=y, xgrid=thetagrid,
                                bwscale=bwscale, smooth=smooth, score=score)
        i <- ccov_np_print_progress(progress=progress, i=i, ii=ii, display=display)
    }

    sirt_progress_cat(progress=progress)
    #-- weights thetagrid
    wgt_thetagrid <- ccov_np_score_density(score=score, thetagrid=thetagrid,
                            smooth=smooth)

    #-- display progress
    if (progress){
        cat('...........................................................\n' )
        cat('Nonparametric Estimation of conditional covariances \n ' )
        utils::flush.console()
    }
    # calculation of conditional covariance
    ccov.table <- data.frame( 'item1ID'=rep( 1:I, I), 'item2ID'=rep( 1:I, each=I ) )
    ccov.table <- ccov.table[ ccov.table$item1ID < ccov.table$item2ID, ]
    ccov.table$N <- apply( ccov.table, 1, FUN=function(ll){
                    sum( rowSums( is.na( data[, c( ll[1], ll[2] ) ] ) )==0 ) } )
    ccov.table <- ccov.table[ ccov.table$N > 0, ]
    ccov.table$item1 <- colnames(data)[ ccov.table$item1ID ]
    ccov.table$item2 <- colnames(data)[ ccov.table$item2ID ]
    ccov.table$itempair <- paste( ccov.table$item1, ccov.table$item2, sep='-' )
    # smoothing all item pairs
    # calculate conditional covariances
    FF <- nrow(ccov.table)
    ccov.matrix <- prod.matrix <- matrix( 0, nrow=length(thetagrid), ncol=FF )
    ii <- 1
    ccov_sum_score <- rep(NA, FF)
    for (ff in 1:FF){
        if (FF>20){
            display <- seq( 1, FF, floor( FF/20 ) )[ 2:20 ]
        } else {
            display <- seq(1,FF)
        }
        data.ff <- data[, c( ccov.table[ff,1], ccov.table[ff,2] ) ]
        which.ff <- which( rowSums( is.na(data.ff) )==0  )
        data.ff <- data.ff[ which.ff, ]
        score.ff <- score[which.ff]
        prod.matrix[,ff] <- ccov_np_regression(x=score.ff, y=data.ff[,1]*data.ff[,2],
                        xgrid=thetagrid, bwscale=bwscale, smooth=smooth, score=score)
        m12 <- icc_items[, ccov.table[ff,1] ]*icc_items[, ccov.table[ff,2] ]
        ccov.matrix[,ff] <- prod.matrix[,ff] - m12

        #- computations based on sum score
        if (use_sum_score){
            res1 <- ccov_np_compute_ccov_sum_score(score=score.ff, data=data.ff)
            score_ff2 <- score.ff - data.ff[,1] - data.ff[,2]
            if (bias_corr){
                res2 <- ccov_np_compute_ccov_sum_score(score=score_ff2, data=data.ff)
            } else {
                res2 <- res1
            }
            ccov_sum_score[ff] <- ( res1$ccov_aggr + res2$ccov_aggr ) / 2
        }

        # print progress
        ii <- ccov_np_print_progress(progress=progress, i=ii, ii=ff, display=display)
    }
    # remove NAs from ccov.matrix
    ccov.matrix[ is.na(ccov.matrix) ] <- 0
    sirt_progress_cat(progress=progress)

    # calculate (weighted) conditional covariance
    ccov.table$ccov <- apply( ccov.matrix, 2, FUN=function(sp){
                                stats::weighted.mean( x=sp, w=wgt_thetagrid ) } )
    if (use_sum_score){
        ccov_sum_score[ is.na(ccov_sum_score) ] <- 0
        ccov.table$ccov <- ccov_sum_score
    }

    #--- output
    res <- list( ccov.table=ccov.table, ccov.matrix=ccov.matrix,
                    data=data, score=score, icc.items=icc_items,
                    wgt.thetagrid=wgt_thetagrid, use_sum_score=use_sum_score,
                    bias_corr=bias_corr, scale_score=scale_score)
    return(res)
}

