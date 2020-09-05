## File Name: expl.detect.R
## File Version: 1.309


#**** Exploratory DETECT analysis
expl.detect <- function( data, score, nclusters, N.est=NULL, seed=NULL,
        bwscale=1.1, smooth=TRUE, use_sum_score=FALSE, hclust_method="ward.D",
        estsample=NULL)
{
    if ( ! is.null(seed) ){
        set.seed(seed)
    }
    # number of items
    I <- ncol(data)
    if (use_sum_score){
        smooth <- FALSE
        score <- rowSums(data)
    }
    scale_score <- TRUE
    if (!smooth){
        scale_score <- FALSE
    }
    # sample for estimation
    N <- nrow(data)
    if ( is.null( N.est ) ){
        N.est <- floor(N/2)
    }
    if (is.null(estsample)){
        estsample <- sort( sample( 1:N, floor( N.est ) ) )
    }
    # validation sample
    valsample <- setdiff( 1:N, estsample )

    #--- Maximizing DETECT index
    # nonparametric estimation of conditional covariance
    ccov_np_args <- list( data=data[ estsample,], score=score[estsample],
                        bwscale=bwscale, smooth=smooth, use_sum_score=use_sum_score,
                        scale_score=scale_score)
    cc_est <- cc <- do.call( what=ccov.np, args=ccov_np_args)
    ccov.matrix <- create.ccov( cc=cc, data=data[ estsample,]  )
    # create distance matrix
    cc1 <- max(ccov.matrix) - ccov.matrix
    # Ward Hierarchical Clustering
    d <- stats::as.dist(cc1)
    fit <- stats::hclust(d, method=hclust_method)         # hierarchical cluster analysis
    clusterfit <- fit
    itemcluster <- data.frame( matrix( 0, I, nclusters ) )
    itemcluster[,1] <- colnames(data)
    colnames(itemcluster) <- c( "item", paste( "cluster", 2:nclusters, sep="") )
    detect.unweighted <- detect.weighted <- NULL
    for (k in 2:nclusters){
        itemcluster[,k] <- stats::cutree( fit, k=k )
        h1 <- detect.index( ccovtable=cc, itemcluster=itemcluster[,k] )
        detect.unweighted <- rbind( detect.unweighted, h1$unweighted )
        detect.weighted <- rbind( detect.weighted, h1$weighted )
    }
    parnames <- c( "DETECT", "ASSI", "RATIO", "MADCOV100", "MCOV100")
    colnames(detect.unweighted) <- paste( parnames, ".est", sep="")
    colnames(detect.weighted) <- paste( parnames, ".est", sep="")
    dfr1 <- data.frame( "N.Cluster"=2:nclusters )
    dfr1$N.items <- I
    dfr1$N.est <- N.est
    dfr1$N.val <- length(valsample)
    dfr1$size.cluster <- sapply( 2:nclusters, FUN=function(tt){
                            paste( table( itemcluster[,tt] ), collapse="-" )
                        } )
    detu <- data.frame( dfr1, detect.unweighted )
    detw <- data.frame( dfr1, detect.weighted )

    #--- Validating DETECT index
    if ( length(valsample) > 0 ){
        ccov_np_args$data=data[ valsample,]
        ccov_np_args$score=score[valsample]
        cc <- do.call( what=ccov.np, args=ccov_np_args)
        detect.unweighted <- detect.weighted <- NULL
        for (k in 2:nclusters){
            h1 <- detect.index( ccovtable=cc, itemcluster=itemcluster[,k] )
            detect.unweighted <- rbind( detect.unweighted, h1$unweighted )
            detect.weighted <- rbind( detect.weighted, h1$weighted )
        }
        colnames(detect.unweighted) <- paste( parnames, ".val", sep="")
        colnames(detect.weighted) <- paste( parnames, ".val", sep="")
        detu <- data.frame( detu, detect.unweighted )
        detw <- data.frame( detw, detect.weighted )
    }
    rownames(detect.unweighted) <- paste0("Cl", 2:nclusters)
    rownames(detect.weighted) <- rownames(detect.unweighted)
    cat("\n\nDETECT (unweighted)\n\n")
    clopt <- which.max( detu$DETECT.est ) + 1
    cat("Optimal Cluster Size is ", clopt, " (Maximum of DETECT Index)\n\n" )
    detu1 <- detu
    for (vv in 6:ncol(detu)){
        detu1[,vv] <- round( detu1[,vv], 3)
    }
    print(detu1)
    res <- list( detect.unweighted=detect.unweighted, detect.weighted=detect.weighted,
                    clusterfit=clusterfit, itemcluster=itemcluster )
    # plot cluster solution
    graphics::plot( res$clusterfit, main=paste( "Cluster Dendogram with ", clopt, " Clusters", sep="") )
    stats::rect.hclust(res$clusterfit, k=clopt, border="red")
    class(res) <- "expl.detect"
    return(res)
}
