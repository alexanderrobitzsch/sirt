## File Name: mcmc_summary.R
## File Version: 0.225



#**** mcmclist descriptives
mcmc_summary <- function( mcmcobj, quantiles=c(.025,.05,.50,.95,.975) )
{
    summary.mcmcobj <- summary(mcmcobj, quantile=quantiles)
    dat.bugs <- mcmc_extract_samples_first_chain(mcmcobj=mcmcobj)
    vars <- colnames(dat.bugs)
    n.iter <- nrow(dat.bugs)
    discret <- round( seq( 1, n.iter, length=4 ) )
    MAP <- Rhat <- rep(0,length(vars) )
    for (ii in seq( 1, length(vars )) ){
        vv <- vars[ii]
        dat.vv <- as.vector( dat.bugs[, vv ] )
        # extract 3 subchains
        l1 <- dat.vv[ discret[1]:discret[2] ]
        l2 <- dat.vv[ (discret[2]+1):discret[3] ]
        l3 <- dat.vv[ (discret[3]+1):discret[4] ]
        W <- ( stats::var(l1) + stats::var(l2) + stats::var(l3) ) / 3
        S1 <- n.iter / 3
        est.chains <- c( mean(l1), mean(l2), mean(l3) )
        B <-  S1 / 2 * sum( ( est.chains - mean( est.chains ) )^2 )
        Rhat[ii] <- ( ( S1-1 ) / S1 * W + B / S1 ) / W
        # mode estimation
        m1 <- stats::density( dat.vv, from=min(dat.vv),
                            to=max(dat.vv) )
        MAP[ii] <- m1$x[ which( m1$y==max( m1$y) ) ]
    }
    res <- data.frame( MAP=MAP, Rhat=Rhat )
    rownames(res) <- vars
    smc3  <- res
    smc2 <- summary.mcmcobj$statistics
    colnames(summary.mcmcobj$quantiles) <- paste0( 'Q', 100*quantiles )
    # calculate effective sample size
    effSize <- sirt_import_coda_effectiveSize( mcmcobj )
    statis <- summary.mcmcobj$statistics
    statis <- cbind( statis[, c(1,2) ],
                apply( as.matrix(mcmcobj), 2, stats::mad ),
                apply( as.matrix(mcmcobj), 2, skewness.sirt ),
                statis[,c(3,4) ]    )
    colnames(statis)[3:4] <- c('MAD', 'skewness' )
    dfr <- data.frame( parameter=rownames(smc3), statis, smc3,
                    SERatio=smc2[,4] / smc2[,2], sampSize=nrow(as.matrix(mcmcobj)),
                    effSize=effSize, summary.mcmcobj$quantiles )
    rownames(dfr) <- NULL
    return(dfr)
}

