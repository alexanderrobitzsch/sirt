## File Name: lsem_weighted_cov.R
## File Version: 0.227
## File Last Change: 2023-03-11

lsem_weighted_cov <- function( x, weights, x_resp=NULL,
        moderator_variable=NULL, loc_linear_smooth=NULL, moderator_value=NULL,
        pd=FALSE)
{
    if (pd){
        requireNamespace('Matrix')
    }
    if (is.null(loc_linear_smooth)){
        loc_linear_smooth <- FALSE
    }
    x <- as.matrix(x)
    if ( is.null(x_resp)){
        x_resp <- 1 - is.na(x)
    }
    eps0 <- 1e-200
    eps <- eps0 * max(weights)
    weights_m <- sqrt( weights + eps ) * x_resp
    x[ ! x_resp ] <- 0

    res <- lsem_weighted_mean( x=x, weights=weights, x_resp=x_resp,
                        moderator_variable=moderator_variable,
                        loc_linear_smooth=loc_linear_smooth,
                        moderator_value=moderator_value)
    x_center <- res$mean
    XC <- matrix( x_center, nrow=nrow(x), ncol=ncol(x), byrow=TRUE )
    x <- x - XC
    weightsN <- crossprod(weights_m)
    xw <- as.matrix( x * weights_m)
    covw <- crossprod(xw) / weightsN
    covw_raw <- covw
    covw2 <- NA*covw

    if (loc_linear_smooth){
        V <- ncol(x)
        for (vv in 1L:V){
            for (ww in vv:V){
                weights1 <- weights
                mod_vv <- stats::lm(x[,vv]~moderator_variable, weights=weights1)
                if (ww==vv){
                    mod_ww <- mod_vv
                } else {
                    mod_ww <- stats::lm(x[,ww]~moderator_variable, weights=weights1)
                }
                rvv <- resid(mod_vv)
                rww <- resid(mod_ww)
                mod <- stats::lm(rvv*rww ~ moderator_variable, weights=weights1)
                cmod <- mod$coefficients
                temp1 <- cmod[1]+cmod[2]*moderator_value
                covw2[vv,ww] <- covw2[ww,vv] <- temp1
            }
        }
        covw <- covw2
    }  # end fit loc_linear_smooth

    if (pd){
        covw <- as.matrix(Matrix::nearPD(x=covw)$mat)
    }

    Nobs <- mean( weightsN[ ! upper.tri(weightsN) ] )

    #-- output
    res <- list( weightsN=weightsN, cov=covw, raw_cov=covw_raw,
                    Nobs=Nobs, mean=x_center)
    return(res)
}
