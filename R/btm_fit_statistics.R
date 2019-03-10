## File Name: btm_fit_statistics.R
## File Version: 0.197


#**** item outfit and infit statistic
btm_fit_statistics <- function( probs, dat0, ind1, ind2, TP, judge=NULL,
    compute_agreement=TRUE)
{
    multiple_judges <- TRUE
    if (is.null(judge)){
        multiple_judges <- FALSE
    }
    # first individual
    X_exp1 <- probs[,1]*1 + 1/2*probs[,3]
    X_var1 <- 1*probs[,1] + 1/4*probs[,3]
    X_var1 <- X_var1 - X_exp1^2
    Z_1 <- ( dat0[,3] - X_exp1 ) / sqrt( X_var1 )

    # second individual
    X_exp2 <- probs[,2]*1 + 1/2*probs[,3]
    X_var2 <- 1*probs[,2] + 1/4*probs[,3]
    X_var2 <- X_var2 - X_exp2^2
    Z_2 <- ( 1 - dat0[,3] - X_exp2 ) / sqrt( X_var2 )

    # compute outfit statistic
    out1 <- rowsum( Z_1^2, dat0[,1] )
    N1 <- rowsum( 1+0*Z_1, dat0[,1] )
    out2 <- rowsum( Z_1^2, dat0[,2] )
    N2 <- rowsum( 1+0*Z_1, dat0[,2] )

    # compute infit statistic
    wvar1 <- rowsum( X_var1, dat0[,1] )
    wvar2 <- rowsum( X_var2, dat0[,2] )
    win1 <- rowsum( X_var1*Z_1^2, dat0[,1] )
    win2 <- rowsum( X_var2*Z_1^2, dat0[,2] )

    out <- btm_fit_combine_tables( win1=out1, win2=out2, ind1=ind1, ind2=ind2, TP=TP )
    N <- btm_fit_combine_tables( win1=N1, win2=N2, ind1=ind1, ind2=ind2, TP=TP )
    wvar <- btm_fit_combine_tables( win1=wvar1, win2=wvar2, ind1=ind1, ind2=ind2, TP=TP )
    win <- btm_fit_combine_tables( win1=win1, win2=win2, ind1=ind1, ind2=ind2, TP=TP )
    outfit <- out/N
    infit <- win/wvar

    #--- fit statistics in case of multiple judges
    fit_judges <- NULL
    if (multiple_judges){
        # outfit and infit statistic
        out1 <- rowsum( Z_1^2, judge )
        N1 <- rowsum( 1+0*Z_1, judge )
        wvar1 <- rowsum( X_var1, judge )
        win1 <- rowsum( X_var1*Z_1^2, judge)
        fit_judges <- data.frame( judge=rownames(out1), outfit=out1/N1,
                            infit=win1/wvar1)

        #* compute agreement statistics
        if (compute_agreement){
            dat1 <- data.frame(judge=judge)
            colnames(dat0) <- c("id1","id2","result")
            ind <- dat0$id1 < dat0$id2
            dat1$id1 <- ifelse(ind, dat0$id1, dat0$id2)
            dat1$id2 <- ifelse(ind, dat0$id2, dat0$id1)
            dat1$result <- ifelse(ind, dat0$result, 1-dat0$result)
            dat1$dyad <- paste0(dat1$id1, "-", dat1$id2)
            a1 <- aggregate( dat1$result, list(dat1$dyad), median )
            colnames(a1) <- c("dyad", "mode")
            a1$N_dyad <- aggregate( 1+0*dat1$result, list(dat1$dyad), sum )[,2]
            a1 <- a1[ ( a1$N_dyad > 2 ) & ( a1$mode %in% c(0,1) ), ]
            dat1$mode <- a1[ match(dat1$dyad, a1$dyad), "mode" ]
            a2 <- aggregate( dat1$result==dat1$mode, list(dat1$judge), mean, na.rm=TRUE)
            fit_judges$agree <- a2[ match(fit_judges$judge, a2[,1]), 2]
        }
    }

    #--- output
    res0 <- list( outfit=outfit, infit=infit, multiple_judges=multiple_judges,
                fit_judges=fit_judges)
    return(res0)
}
