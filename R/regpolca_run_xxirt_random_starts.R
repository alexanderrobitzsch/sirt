## File Name: regpolca_run_xxirt_random_starts.R
## File Version: 0.14

regpolca_run_xxirt_random_starts <- function(args, random_starts, sd_noise_init,
    par_item_init=NULL)
{
    partable1 <- partable <- args$partable
    customTheta1 <- customTheta <- args$customTheta
    NR <- nrow(partable)
    CP <- length(customTheta$par)
    opt0 <- Inf
    partable_opt <- NULL
    par_Theta_opt <- NULL
    rr_opt <- NULL
    if (!is.null(par_item_init)){
        args$partable$value <- par_item_init
    }
    if (random_starts>0){
        for (rr in 1L:random_starts){
            fac_rr <- sqrt(rr-1)
            sd_rr <- sd_noise_init*fac_rr/random_starts
            partable1$value <- partable$value + partable$est*stats::rnorm(NR, sd=sd_rr)
            args$verbose_index <- paste0("Random start ",rr)
            args$partable <- partable1
            customTheta1$par <- customTheta$par + stats::rnorm( CP, sd=sd_rr )
            args$customTheta <- customTheta1
            cat(paste0( "*** Random start ", rr, " ***"), "\n")
            mod <- do.call(what=xxirt, args=args)
            par_Theta0 <- mod$customTheta$par
            partable0 <- mod$partable
            opt_val <- mod$opt_val
            if (opt_val < opt0){
                opt0 <- opt_val
                par_Theta_opt <- par_Theta0
                partable_opt <- partable0
                rr_opt <- rr
            }
        }
        partable1$value <- partable_opt$value
        args$partable <- partable1
        customTheta1$par <- par_Theta_opt
        args$customTheta <- customTheta1
    }
    return(args)
}
