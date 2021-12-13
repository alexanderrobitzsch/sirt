## File Name: linking_haberman_als.R
## File Version: 0.654



#--- alternating least squares for Haberman linking
linking_haberman_als <- function(logaM, wgtM, maxiter, conv,
            progress, est.type, cutoff, reference_value=0, adjust_main_effects=FALSE,
            estimation="OLS", lts_prop=.5, vcov=TRUE)
{
    non_null <- sum( abs(logaM) > 0, na.rm=TRUE)
    iter <- 0
    parchange <- 1000
    NS <- ncol(logaM)
    NI <- nrow(logaM)
    #-- initial values study parameters
    logaAt <- rep(0,NS)
    At_inits <- TRUE
    if ( At_inits ){
        logaAt <- weighted_colMeans( mat=logaM, wgt=wgtM )
        logaAt <- logaAt - logaAt[1]
    }
    wgtM0 <- wgtM
    wgt_adj <- 1 + 0*wgtM
    eps <- 1E-5
    wgtM <- wgtM + eps
    #- initial OLS estimation of item parameters
    logaAt_M <- sirt_matrix2( x=logaAt, nrow=NI)
    logaM_adj1 <- logaM - logaAt_M
    logaj <- weighted_rowMeans( mat=logaM_adj1, wgt=wgtM )
    if (! ( estimation %in% c("L0","BSQ","HUB") ) ) {
        cutoff <- Inf
    }
    if (estimation=="L0"){
        res1 <- L0_polish(x=logaM, tol=cutoff, type=1)
        wgtM0 <- wgtM <- res1$wgt
        logaM1 <- res1$x_update
    }
    I <- nrow(wgtM)
    NS <- ncol(wgtM)
    abs_a_change <- 1

    #*** begin algorithm
    while( ( parchange > conv ) & (iter < maxiter) ){
        logaAt0 <- logaAt
        logaj_old <- logaj
        #--- estimate item parameters
        res <- linking_haberman_als_residual_weights( logaj=logaj, logaAt=logaAt,
                    logaM=logaM, cutoff=cutoff, wgtM0=wgtM0, eps=eps,
                    estimation=estimation, lts_prop=lts_prop)
        loga_resid <- res$loga_resid
        wgtM <- res$wgtM
        wgt_adj <- res$wgt_adj

        # estimation
        logaAt_M <- sirt_matrix2( x=logaAt, nrow=NI)
        logaM_adj1 <- logaM - logaAt_M
        logaj <- rep(NA,I)
        if (estimation %in% c("OLS","BSQ","HUB")){
            logaj <- weighted_rowMeans( mat=logaM_adj1, wgt=wgtM )
        }
        if (estimation %in% c("MED")){
            for (ii in 1:I){
                logaj[ii] <- linking_haberman_compute_median(x=logaM_adj1[ii,], w=wgtM[ii,])
            }
        }
        if (estimation %in% c("L0","L1")){
            if (estimation=="L1"){
                logaM1 <- logaM
            }
            res1 <- L1_polish(x=logaM1, type=1)
            logaj <- res1$row
            logaAt <- res1$col
        }
        if (estimation %in% c("LTS")){
            for (ii in 1:I){
                logaj[ii] <- linking_haberman_compute_lts_mean(x=logaM_adj1[ii,], w=wgtM[ii,],
                                    lts_prop=lts_prop)
            }
        }

        #--- estimate linking parameters
        logaMadj <- logaM - logaj
        res <- linking_haberman_als_residual_weights( logaj=logaj, logaAt=logaAt,
                    logaM=logaM, cutoff=cutoff, wgtM0=wgtM0, eps=eps,
                    estimation=estimation, lts_prop=lts_prop)
        wgtM <- res$wgtM
        wgt_adj <- res$wgt_adj
        loga_resid <- res$loga_resid
        cutoff_used <- res$cutoff
        k_estimate <- res$k_estimate

        #* estimation of parameters
        if (estimation %in% c("OLS","BSQ","LTS","HUB")){
            logaAt <- weighted_colMeans( mat=logaMadj, wgt=wgtM )
        }
        if (estimation %in% c("MED")){
            for (ss in 1:NS){
                logaAt[ss] <- linking_haberman_compute_median(x=logaMadj[,ss], w=wgtM[,ss])
            }
        }
        if ( ! adjust_main_effects){
            logaAt[1] <- reference_value
        } else {
            ma <- logaAt[1]
            logaAt <- logaAt - ma + reference_value
        }
        a_change <- logaAt - logaAt0

        #- stabilize convergence
        if (iter>50){
            if (max(abs(a_change)) >=abs_a_change ){
                a_change <- ifelse( abs(a_change) >=abs_a_change, .95*abs_a_change, a_change )
                logaAt <- logaAt0 + a_change
            }
        }
        parchange <- abs_a_change <- max(abs(a_change))
        if (progress){
            cat( paste0( "** ", est.type, " estimation | Iteration ", iter+1, " | ",
                "Max. parameter change=", round( parchange, 6 ) ), "\n")
            utils::flush.console()
        }
        iter <- iter + 1

        #-- stop iterations
        if (estimation %in% c("L0","L1")){
            break
        }

    }
    if (progress){
        cat("\n")
    }

    #------- summary of regression

    # residual SD
    selitems <- which( rowSums( 1 - is.na( loga_resid ) ) > 1 )

    #--- calculation of standard errors of regression coefficients
    sd_loga <- stats::sd(logaAt)
    if ( sd_loga < 1E-10 ){
        res <- list( vcov=0*diag(NS-1), se=rep(0,NS-1)  )
    } else {
        res <- linking_haberman_als_vcov( regr_resid=loga_resid,
                    regr_wgt=wgtM, selitems=selitems, transf_pars=logaAt,
                    estimation=estimation, vcov=vcov, NS=NS)
    }
    #--- item statistics
    item_stat <- data.frame( study=colnames(wgtM0) )
    item_stat$N_items <- colSums( wgtM0 > 0, na.rm=TRUE)
    item_stat$sumwgt_items <- colSums( (wgtM0 > 0)*wgt_adj, na.rm=TRUE )
    if (estimation!="LTS"){
        lts_prop <- 1
    }
    #*** end algorithm
    res <- list( logaAt=logaAt, logaj=logaj, loga_resid=loga_resid, loga_wgt=wgtM,
                loga_wgt_adj=wgt_adj, vcov=res$vcov, se=c(NA, res$se),
                item_stat=item_stat, iter=iter, cutoff=cutoff_used, estimation=estimation,
                k_estimate=k_estimate, lts_prop=lts_prop)
    return(res)
}

