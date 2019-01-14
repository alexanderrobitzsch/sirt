## File Name: linking.haberman.R
## File Version: 2.635


#**** Linking Haberman: ETS Research Report 2009
linking.haberman <- function( itempars, personpars=NULL,
        estimation="OLS", a_trim=Inf, b_trim=Inf, lts_prop=.5,
        a_log=TRUE, conv=.00001, maxiter=1000, progress=TRUE, adjust_main_effects=TRUE)
{
    CALL <- match.call()
    s1 <- Sys.time()

    #--- process item parameters
    res <- linking_proc_itempars(itempars=itempars)
    itempars <- res$itempars
    NS <- res$NS
    NI <- res$NI
    items <- res$items
    studies <- res$studies
    wgtM <- res$wgtM
    aM <- res$aM
    bM <- res$bM
    est_pars <- res$est_pars
    weights_exist <- res$weights_exist

    a.orig <- aM
    b.orig <- bM

    #**** estimation of A
    if (a_log){
        logaM <- log(aM)
    } else {
        logaM <- aM
    }
    est_type <- "A (slopes)"
    resA <- linking_haberman_als(logaM=logaM, wgtM=wgtM, maxiter=maxiter, conv=conv,
                    progress=progress, est.type=est_type, cutoff=a_trim,
                    reference_value=1-a_log, adjust_main_effects=adjust_main_effects,
                    estimation=estimation, lts_prop=lts_prop)
    if (a_log){
        aj <- exp(resA$logaj)
        At <- exp(resA$logaAt)
    } else {
        aj <- resA$logaj
        At <- resA$logaAt
    }
    aj_resid <- resA$loga_resid
    aj_wgt_adj <- resA$loga_wgt_adj
    aj_wgtM <- resA$loga_wgt
    aj_vcov <- resA$vcov
    aj_se <- resA$se
    aj_item_stat <- resA$item_stat

    ## vcov of transformed parameters
    # H1 <- diag( At[-1] )
    H1 <- diag2( At[-1] )
    res <- linking_haberman_vcov_transformation( H1=H1, aj_vcov=aj_vcov )
    aj_vcov <- res$vcov
    aj_se <- c( NA, res$se )

    #**** estimation of B
    est_type <- "B (intercepts)"
    At_m <- sirt_matrix2( At, nrow=NI)
    bMadj <- bM * At_m
    resB <- linking_haberman_als(logaM=bMadj, wgtM=wgtM, maxiter=maxiter,
                    conv=conv, progress=progress, est.type=est_type,
                    cutoff=b_trim, reference_value=0, adjust_main_effects=adjust_main_effects,
                    estimation=estimation, lts_prop=lts_prop)
    Bj <- resB$logaj
    Bt <- resB$logaAt
    Bj_resid <- resB$loga_resid
    Bj_wgt_adj <- resB$loga_wgt_adj
    Bj_wgtM <- resB$loga_wgt
    Bj_vcov <- resB$vcov
    Bj_se <- resB$se
    Bj_item_stat <- resB$item_stat

    #*****
    # transformations
    transf.pars <- data.frame( study=studies, At=At, se_At=aj_se, Bt=Bt, se_Bt=Bj_se )
    rownames(transf.pars) <- NULL
    transf.itempars <- data.frame( study=studies, At=1/At, se_At=NA, At2=At,
            se_At2=transf.pars$se_At, Bt=Bt, se_Bt=transf.pars$se_Bt )
    rownames(transf.itempars) <- NULL

    #*** transformation 1/At
    H1 <- diag2( - 1 / ( At[-1] )^2 )

    res <- linking_haberman_vcov_transformation( H1=H1, aj_vcov=aj_vcov )
    transf.itempars$se_At <- c( NA, res$se )

    #**** tranbsformation for person parameters
    transf.personpars <- transf.itempars[,c("study","At","se_At2","Bt", "se_Bt")]
    transf.personpars$At <- transf.pars$At
    transf.personpars$Bt <- -transf.pars$Bt
    colnames(transf.personpars) <- c("study", "A_theta",
                "se_A_theta", "B_theta", "se_B_theta" )
    colnames(transf.itempars) <- c("study", "A_a", "se_A_a",
                "A_b", "se_A_b", "B_b", "se_B_b" )
    # new item parameters
    joint.itempars <- data.frame("item"=items, "aj"=aj, "bj"=Bj )
    # transformations for item parameters
    AtM <- sirt_matrix2( At, nrow=NI)
    BtM <- sirt_matrix2( Bt, nrow=NI)
    aM <- aM / AtM
    bM <- bM * AtM - BtM
    #****
    # transform person parameters
    if ( ! is.null( personpars) ){
        for (ll in 1:NS){
            pp0 <- pp1 <- personpars[[ll]]
            pp1 <- transf.personpars$A_theta[ll] * pp1 + transf.personpars$B_theta[ll]
            ind <- which( substring( colnames(pp0),1,2) %in% c("se", "SE") )
            if ( length(ind) > 0 ){
                pp1[,ind] <- transf.personpars$A_theta[ll] * pp0[,ind]
            }
            ind <- which( substring( colnames(pp0),1,3) %in% c("pid") )
            if ( length(ind) > 0 ){
                pp1[,ind] <- pp0[,ind]
            }
            personpars[[ll]] <- pp1
        }
    }

    #****** calculate R-squared measures of invariance
    # select items for R2 calculation for which at least
    # two studies are available.
    selitems <- which( rowSums( 1 - is.na( a.orig ) ) > 1 )
    Rsquared.partial.invariance <- Rsquared.invariance <- c(NA,NA)
    names(Rsquared.invariance) <- c("slopes", "intercepts" )
    names(Rsquared.partial.invariance) <- c("slopes", "intercepts" )

    # retransformed parameters
    aj1 <- aj * AtM
    a.res <- a.orig - aj1

    Rsquared.invariance["slopes"] <- 1 - sirt_sum( a.res[ selitems,]^2 ) / sirt_sum( a.orig[ selitems, ]^2 )
    Rsquared.partial.invariance["slopes"] <- 1 -
        sirt_sum( a.res[ selitems,]^2 * aj_wgtM[selitems,]  ) /
        sirt_sum( a.orig[ selitems, ]^2 * aj_wgtM[selitems, ]  )
    bj1 <- 1 / AtM *( Bj + BtM )
    b.res <- b.orig - bj1

    Rsquared.partial.invariance["intercepts"] <- 1 -
        sirt_sum( b.res[ selitems,]^2 * Bj_wgtM[selitems,] ) /
        sirt_sum( b.orig[ selitems, ]^2 * Bj_wgtM[selitems,] )
    Rsquared.invariance["intercepts"] <- 1 -
        sirt_sum( b.res[ selitems,]^2  ) / sirt_sum( b.orig[ selitems, ]^2 )
    es.invariance <- rbind( Rsquared.invariance,
            sqrt( 1- Rsquared.invariance ) )
    rownames(es.invariance) <- c("R2", "sqrtU2")
    es.partial.invariance <- rbind( Rsquared.partial.invariance,
            sqrt( 1- Rsquared.partial.invariance ) )
    rownames(es.partial.invariance) <- c("R2", "sqrtU2")

    # logical indicating whether item slopes are linked
    linking_slopes <- stats::sd( transf.pars$At ) < 1E-10

    s2 <- Sys.time()
    time <- list(s1=s1, s2=s2)

    #--- output
    res <- list( transf.pars=transf.pars, transf.itempars=transf.itempars,
            transf.personpars=transf.personpars, joint.itempars=joint.itempars,
            a.trans=aM, b.trans=bM, a.orig=a.orig, b.orig=b.orig,
            a.resid=aj_resid, b.resid=Bj_resid, personpars=personpars,
            es.invariance=es.invariance, es.robust=es.partial.invariance,
            selitems=selitems, a_trim=a_trim, b_trim=b_trim,
            a.wgt=aj_wgtM, b.wgt=Bj_wgtM, a.wgt.adj=aj_wgt_adj, b.wgt.adj=Bj_wgt_adj,
            a.vcov=aj_vcov, b.vcov=Bj_vcov, a.item_stat=aj_item_stat,
            b.item_stat=Bj_item_stat, linking_slopes=linking_slopes,
            description='Linking according to Haberman (2009)',
            res_opt_slopes=resA, res_opt_intercepts=resB, CALL=CALL, time=time )
    class(res) <- "linking.haberman"
    return(res)
}


