## File Name: lsem_weighted_mean.R
## File Version: 0.173

lsem_weighted_mean <- function( x, weights, x_resp=NULL,
        moderator_variable=NULL, loc_linear_smooth=NULL, moderator_value=NULL )
{
    if (is.null(loc_linear_smooth)){
        loc_linear_smooth <- FALSE
    }
    x <- as.matrix(x)
    if ( is.null(x_resp)){
        x_resp <- 1 - is.na(x)
    }
    weights_m <- weights * x_resp
    x[ ! x_resp ] <- 0
    weightsN <- colSums(weights_m)
    wm_raw <- wm <- colSums( x * weights_m ) / weightsN

    if (loc_linear_smooth){
        V <- length(wm)
        wm2 <- rep(NA, V)
        colnames(wm2) <- colnames(wm)
        for (vv in 1L:V){
            mod <- stats::lm(x[,vv]~moderator_variable, weights=weights_m[,vv])
            cmod <- mod$coefficients
            wm2[vv] <- cmod[1]+cmod[2]*moderator_value
        }
        wm <- wm2
    }

    #-- output
    res <- list( weightsN=weightsN, mean=wm, raw_mean=wm_raw )
    return(res)
}
