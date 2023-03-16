## File Name: mcmc_WaldTest.R
## File Version: 0.25


#** Wald Test for a set of hypotheses
mcmc_WaldTest <- function( mcmcobj, hypotheses )
{
    mcmcobj <- mcmc_extract_samples_first_chain(mcmcobj=mcmcobj)
    NH <- length(hypotheses)
    n1 <- ncol(mcmcobj)
    mcmcobj <- mcmc_derivedPars( mcmcobj=mcmcobj, derivedPars=hypotheses)
    n2 <- ncol(mcmcobj)
    mcmcobj <- mcmcobj[, seq(n1+1,n2) ]
    v1 <- mcmc_vcov(mcmcobj)
    s1 <- mcmc_summary(mcmcobj)
    c1 <- s1$MAP
    # compute test statistic
    v1a <- sirt_import_MASS_ginv(X=v1)
    v1_svd <- svd(v1)
    eps <- 1E-10
    NH <- sum( v1_svd$d > eps )
    W <- t(c1) %*% v1a %*% c1
    stat <- c( "chi2"=W, "df"=NH)
    stat["p"] <- 1 - stats::pchisq( W, df=stat["df"])
    res <- list( hypotheses_summary=s1, chisq_stat=stat)
    class(res) <- "mcmc_WaldTest"
    return(res)
}

