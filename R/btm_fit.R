## File Name: btm_fit.R
## File Version: 0.13


############################################
# item outfit and infit statistic
btm_fit <- function( probs, dat0, ind1, ind2, TP)
{
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
    wvar1 <- rowsum( X_var1, dat0[,1] )
    wvar2 <- rowsum( X_var2, dat0[,2] )
    win1 <- rowsum( X_var1*Z_1^2, dat0[,1] )
    win2 <- rowsum( X_var2*Z_1^2, dat0[,2] )

    out <- btm_fit_combine_tables( win1=out1, win2=out2, ind1=ind1, ind2=ind2, TP=TP )
    N <- btm_fit_combine_tables( win1=N1, win2=N2, ind1=ind1, ind2=ind2, TP=TP )
    wvar <- btm_fit_combine_tables( win1=wvar1, win2=wvar2, ind1=ind1, ind2=ind2, TP=TP )
    win <- btm_fit_combine_tables( win1=win1, win2=win2, ind1=ind1, ind2=ind2, TP=TP )
    #--- output
    res0 <- list( outfit=out / N,  infit=win / wvar )
    return(res0)
}
#########################################################################

