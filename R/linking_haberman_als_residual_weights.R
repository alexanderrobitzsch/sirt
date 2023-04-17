## File Name: linking_haberman_als_residual_weights.R
## File Version: 0.372


linking_haberman_als_residual_weights <- function( logaj, logaAt,
        logaM, cutoff, wgtM0, eps, estimation="OLS", lts_prop=.5 )
{
    loga_expected <- TAM::tam_outer( x=logaj, y=logaAt, op='+' )
    loga_resid <- logaM - loga_expected
    NS <- ncol(wgtM0)
    NI <- nrow(wgtM0)
    wgt_adj <- matrix(1, nrow=NI, ncol=NS)
    wgtM <- wgtM0
    k_estimate <- FALSE
    eps_k <- 1e-5

    #*--- settings
    if (estimation=='BSQ'){
        k_fac <- 4.685
        wgt_fn <- linking_haberman_bisquare_weight
    }
    if (estimation=='HUB'){
        k_fac <- 1.345
        wgt_fn <- linking_haberman_huber_weight
    }

    #********************
    #-- estimation BSQ or HUB
    if (estimation %in% c('BSQ','HUB') ){
        if (cutoff==Inf ){
            k <- k_fac*mad_normalized(x=loga_resid)
            k_estimate <- TRUE
            cutoff <- k + eps_k
        }
        if (cutoff<Inf){
            min_x <- apply(abs(loga_resid), 2, min, na.rm=TRUE)
            min_x <- 2*max(min_x)
            cutoff <- max(cutoff, min_x)
            args <- list(x=loga_resid, cutoff=cutoff)
            wgt_adj <- do.call(wgt_fn, args=args )
        }
    }
    #-- estimation LTS
    if (estimation=='LTS'){
        for (ss in 1:NS){
            e <- loga_resid[,ss]
            e <- e - median(e, na.rm=TRUE)
            dfr_resid <- data.frame(item=1:NI, e=e )
            dfr_resid <- na.omit(dfr_resid)
            dfr_resid <- dfr_resid[ order(abs(dfr_resid$e), decreasing=TRUE), ]
            wgt_adj[ is.na(loga_resid[,ss]), ss ] <- 0
            n <- nrow(dfr_resid)
            n_del <- floor( (1-lts_prop)*n)
            m1 <- dfr_resid[ 1:n_del, c(2,1) ]
            wgt_adj[ dfr_resid[ 1:n_del, 1 ], ss ] <- 0
        }
    }

    #*** update weights
    wgt_adj[ is.na(wgt_adj) ] <- 0
    wgtM <- wgtM0 * wgt_adj + eps

    #--- output
    res <- list(loga_resid=loga_resid, wgt_adj=wgt_adj, wgtM=wgtM,
                    cutoff=cutoff, k_estimate=k_estimate)
    return(res)
}

