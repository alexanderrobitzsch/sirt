## File Name: noharm_sirt_compute_chi_square_statistics.R
## File Version: 0.04

noharm_sirt_compute_chi_square_statistics <- function(res, residuals, pm,
    N, I, sumwgtm, modesttype, Nestpars)
{
    if (modesttype==1){
        RM <- residuals
        PV <- diag(pm)
        g1 <- sqrt( TAM::tam_outer( PV * (1-PV), PV * ( 1-PV ) ) )
        rM <- ( RM / g1  )
        zM <- 0.5 * log( 1 + rM ) - 0.5 * log( 1 - rM )
        # chi square
        res$chisquare <- X2 <- ( N - 3 ) * sum( zM^2 )
        res$df <- df <- I + sumwgtm - Nestpars$total
        res$chisquare_df  <- res$chisquare / res$df
        # calculate RMSEA
        res$rmsea <- rmsea <- sqrt(max(c( (X2 / res$Nobs ) / df - 1/res$Nobs, 0)))
        # calculate p values
        res$p.chisquare <- 1 - stats::pchisq( res$chisquare, df=res$df )
    }
    #--- output
    return(res)
}
