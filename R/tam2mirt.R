## File Name: tam2mirt.R
## File Version: 0.292
## File Last Change: 2023-03-08


# convert a fitted tam object into a mirt object
tam2mirt <- function( tamobj )
{
    est.mirt <- FALSE
    # extract intercept
    AXsi <- tamobj$AXsi
    # extract loadings
    B <- tamobj$B
    # number of dimensions
    D <- dim(B)[3]
    # extract trait distribution
    mean.trait <- tamobj$beta
    cov.trait <- tamobj$variance
    # extract data
    dat <- tamobj$resp
    # factors
    if (D==1){
        factors <- 'F1'
    }
    if (D>1){
        factors <- dimnames(tamobj$B)[[3]]
    }
    # lavaan syntax with fixed values
    lavsyn <- tam2mirt_fix( D=D, factors=factors, B=B, dat=dat, AXsi=AXsi,
                    mean.trait=mean.trait, cov.trait=cov.trait, tamobj=tamobj )
    # lavaan syntax with freed values
    lavsyn.freed <- tam2mirt_freed( D=D, factors=factors, B=B, dat=dat,
            AXsi=AXsi, mean.trait=mean.trait, cov.trait=cov.trait, tamobj=tamobj )
    # pseudo-estimate model in mirt: just create mirt object structure
    res <- lavaan2mirt( dat=dat, lavmodel=lavsyn, est.mirt=TRUE )
    #--- include parameters in mirt object
    res$mirt@Model$nest <- as.integer(tamobj$ic$np ) # number of estimated parameters
    # recalculate AIC, BIC, AICc and SABIC
    res$mirt@Fit$AIC <- tamobj$ic$AIC
    res$mirt@Fit$BIC <- tamobj$ic$BIC
    res$mirt@Fit$AICc <- tamobj$ic$AICc
    res$mirt@Fit$SABIC <- tamobj$ic$aBIC
    # use theta grid from estimation in TAM
    res$mirt@Model$Theta <- tamobj$theta
    res$mirt@Options$quadpts <- nrow(tamobj$theta)
    # output
    res$lavaan.syntax.fixed <- lavsyn
    res$lavaan.syntax.freed <- lavsyn.freed
    # res$tamobj <- tamobj
    return(res)
}
