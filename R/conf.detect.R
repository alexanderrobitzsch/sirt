## File Name: conf.detect.R
## File Version: 1.211


# Confirmatory DETECT analysis
conf.detect <- function( data, score, itemcluster, bwscale=1.1, progress=TRUE,
                thetagrid=seq( -3,3,len=200), smooth=TRUE, use_sum_score=FALSE,
                bias_corr=TRUE)
{
    CALL <- match.call()
    cat('-----------------------------------------------------------\n' )
    cat('Confirmatory DETECT Analysis \n' )
    utils::flush.console()
    h1 <- is.matrix(score)
    if (h1){
        PP <- ncol(score)
    }
    is_one_score <- TRUE
    if (! h1 ){
        cat('Conditioning on 1 Score\n' )
    } else {
        cat(paste('Conditioning on ',PP, ' Scores\n', sep='') )
        is_one_score <- FALSE
    }
    cat(paste('Bandwidth Scale:', bwscale, '\n' ) )
    utils::flush.console()
    scale_score <- TRUE
    if (!smooth){
        scale_score <- FALSE
    }
    args_ccov_np <- list( data=data, score=score, bwscale=bwscale,
                            progress=progress, thetagrid=thetagrid, smooth=smooth,
                            scale_score=scale_score, use_sum_score=use_sum_score,
                            bias_corr=bias_corr)
    if ( ! h1 ){
        ccovtable <- do.call( what=ccov.np, args=args_ccov_np)
        res <- detect.index( ccovtable=ccovtable, itemcluster=itemcluster )
    } else {
        ccovtable.list <- list()
        args_ccov_np$progress <- FALSE
        for (pp in 1:PP){
            cat( paste( 'DETECT Calculation Score ', pp, '\n', sep='') ) ;
            utils::flush.console()
            args_ccov_np$score <- score[,pp]
            ccovtable.list[[pp]] <- do.call( what=ccov.np, args=args_ccov_np)
        }
        detect.list <- lapply( ccovtable.list, FUN=function(ccovtable){
                    detect.index( ccovtable, itemcluster=itemcluster ) } )
        detect.matrix <- matrix( unlist( lapply( detect.list, FUN=function( ll){
                                c( ll[1,], ll[2,], ll[3,] ) } ) ), nrow=PP, byrow=TRUE)
        detect.summary <- data.frame( NScores=PP, Mean=colMeans( detect.matrix ),
                                SD=apply( detect.matrix, 2, stats::sd ),
                                Min=apply( detect.matrix, 2, min ),
                                Max=apply( detect.matrix, 2, max ) )
        rownames(detect.summary) <- c('DETECT Unweighted', 'DETECT Weighted',
                                        'ASSI Unweighted', 'ASSI Weighted',
                                        'RATIO Unweighted', 'RATIO Weighted' )
    }
    cat('-----------------------------------------------------------\n' )
    if ( ! h1){
        res <- list( detect=res, ccovtable=ccovtable, detect.summary=res )
    } else {
        res <- list( detect=detect.list, ccovtable=ccovtable.list,
                        detect.summary=detect.summary )
    }
    res$is_one_score <- is_one_score
    res$CALL <- CALL
    res$bwscale <- bwscale
    res$itemcluster <- itemcluster
    res$smooth <- smooth
    #--- print
    print(round(res$detect.summary,3))
    #--- return
    class(res) <- 'conf.detect'
    return(res)
}

