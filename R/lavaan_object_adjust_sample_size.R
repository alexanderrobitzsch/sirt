## File Name: lavaan_object_adjust_sample_size.R
## File Version: 0.106

#- adjust standard errors in lavaan object for modified sample size
lavaan_object_adjust_sample_size <- function(object, n_used)
{
    object2 <- object

    # adjustment factor for standard errors
    Data <- object@Data
    ngroups <- Data@ngroups
    nobs_all <- sum(unlist(Data@nobs))
    fac2 <- nobs_all / n_used
    fac <- sqrt(fac2)

    #- Data@nobs
    nobs <- nobs0 <- Data@nobs
    for (gg in 1:ngroups){
        nobs[[gg]] <- round(nobs[[gg]] / fac2, 0)
    }
    Data@nobs <- nobs

    #- Data@weights
    weights <- Data@weights
    NW <- unlist(weights)
    if (is.null(NW)){
        weights <- list()
        for (gg in 1:ngroups){
            weights[[gg]] <- rep(1, nobs0)
        }
    }
    for (gg in 1:ngroups){
        weights[[gg]] <- weights[[gg]] / fac2
    }
    Data@weights <- weights

    #- ParTable$se
    ParTable <- object@ParTable
    ParTable$se <- ParTable$se * fac

    #- Fit$se
    Fit <- object@Fit
    Fit@se <- Fit@se * fac

    #- vcov@vcov
    vcov <- object@vcov
    vcov$vcov <- vcov$vcov * fac2

    #** include all redefined objects
    object2@Data <- Data
    object2@ParTable <- ParTable
    object2@Fit <- Fit
    object2@vcov <- vcov

    #-- output
    return(object2)
}
