## File Name: tam2mirt_fix.R
## File Version: 0.14



##################################################################
# return lavaan syntax with fixed parameters
tam2mirt_fix <- function( D, factors, B, dat, AXsi,
        mean.trait, cov.trait, tamobj )
{
    # create lavaan syntax with constraints
    lavsyn <- NULL
    for (dd in 1:D){
        # dd <- 1
        fac.dd <- factors[dd]
        # create terms for loadings
        B2.dd <- round( B[,2,dd], 4)
        syn0 <- paste0( paste0( B2.dd[ B2.dd!=0], "*", colnames(dat)[ B2.dd!=0] ), collapse="+" )
        syn0 <- paste0( fac.dd, "=~ ", syn0, "\n")
        lavsyn <- paste0( lavsyn, syn0 )
    }
    # create syntax for intercepts
    maxK <- ncol(AXsi) - 1
    for (kk in 1:maxK){
        t1 <- round( AXsi[,kk+1], 4 )
        string1 <- paste0("t", kk )
        syn0 <- paste0( colnames(dat), " | ", t1, "*", string1)
        syn0 <- paste0( syn0, collapse="\n")
        hh <- ""
        if (kk !=maxK){ hh <- "\n" }
        lavsyn <- paste0( lavsyn, syn0, hh)
    }
    # guessing and slipping parameters
    itemg <- colnames(dat)[ maxK==1 ]
    if ( length(itemg) > 0 ){
        lavsyn <- paste0( lavsyn, "\n", paste0( paste0( itemg, " ?=0*g1" ), collapse="\n") )
        lavsyn <- paste0( lavsyn, "\n", paste0( paste0( itemg, " ?=0*s1" ), collapse="\n") )
    }
    # syntax for means
    syn0 <- paste0( factors, " ~ ", round(as.vector(mean.trait),4), "*1"  )
    syn0 <- paste0( syn0, collapse="\n")
    lavsyn <- paste0( lavsyn, "\n", syn0 )
    # syntax for variances
    syn0 <- paste0( factors, " ~~ ", round( as.vector(diag(cov.trait)),4), "*",factors  )
    syn0 <- paste0( syn0, collapse="\n")
    lavsyn <- paste0( lavsyn, "\n", syn0 )
    # syntax for covariances
    if (D>1){
        for (dd in 1:(D-1)){
            for (ee in (dd+1):D ){
                syn0 <- paste0( factors[dd], " ~~ ",
                            round( cov.trait[dd,ee],4), "*",factors[ee]  )
                syn0 <- paste0( syn0, collapse="\n")
                lavsyn <- paste0( lavsyn, "\n", syn0 )
            }
        }
    }
    # finalize lavaan syntax
    lavsyn <- paste0( lavsyn, " \n")
    return(lavsyn)
}
##################################################################
