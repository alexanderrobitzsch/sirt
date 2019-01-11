## File Name: linking_haberman_als_residual_weights.R
## File Version: 0.362


linking_haberman_als_residual_weights <- function( logaj, logaAt,
        logaM, cutoff, wgtM0, eps, estimation="OLS", lts_prop=.5 )
{
    loga_expected <- TAM::tam_outer( x=logaj, y=logaAt, op="+" )
    loga_resid <- logaM - loga_expected
    NS <- ncol(wgtM0)    
    NI <- nrow(wgtM0)
    wgt_adj <- matrix(1, nrow=NI, ncol=NS)    
    wgtM <- wgtM0
    k_estimate <- FALSE
    
    #-- estimation BSQ
    if (estimation == "BSQ"){
        if (cutoff==Inf ){
            eps_k <- 1e-5
            k <- linking_haberman_bisquare_regression_tuning_constant(x=loga_resid)
            k_estimate <- TRUE
            cutoff <- k + eps_k
        }
        if (cutoff<Inf){            
            min_x <- apply(abs(loga_resid), 2, min, na.rm=TRUE)        
            min_x <- 2*max(min_x)
            cutoff <- max(cutoff, min_x)
            wgt_adj <- ( 1 - ( loga_resid / cutoff )^2 )^2
            wgt_adj <- (abs(loga_resid) <=cutoff)*wgt_adj            
        }
    }
    #-- estimation LTS
    if (estimation=="LTS"){        
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

