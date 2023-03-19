## File Name: lsem.test.R
## File Version: 0.133

#**** test LSEM model based on bootstrap
lsem.test <- function( mod, bmod, models=NULL )
{
    parameters <- mod$parameters
    repl_factor <- bmod$repl_factor
    R <- bmod$R
    if (is.null(repl_factor)){
        repl_factor <- 1/(R-1)
    }

    parnames <- unique(paste(parameters$par))
    w <- mod$moderator.density$wgt
    NG <- length(w)
    bmod_missing <- missing(bmod)
    if (bmod_missing){
        parameters_boot <- NULL
    } else {
        parameters_boot <- bmod$parameters_boot
    }

    #* define design matrix
    A <- matrix(0, nrow=NG-1, ncol=NG)
    for (gg in 1:(NG-1)){
        A[gg,gg] <- -1
        A[gg,gg+1] <- 1
    }

    wald_test_global <- NULL
    dfr <- NULL

    if (! bmod_missing){
        for (pp in parnames){
            ind_pp <- which( parameters$par==pp )
            parameters_pp <- parameters[ind_pp, ]
            theta <- parameters_pp$est
            par_boot_pp <- t( parameters_boot[ind_pp, ] )
            V <- stats::cov( par_boot_pp )*(R-1)*repl_factor
            M <- TAM::weighted_mean(x=theta, w=w)
            SD <- TAM::weighted_sd(x=theta, w=w)
            dfr1 <- data.frame(par=pp, M=M, SD=SD, chisq=NA, df=NA, p=NA)
            if (SD>1e-10){
                res <- lsem_wald_test(theta=theta, V=V, A=A)
                dfr1$chisq <- res$chisq
                dfr1$df <- res$df
                dfr1$p <- res$p
            }
            dfr <- rbind(dfr, dfr1)
        }
        wald_test_global <- dfr
    }


    #** run models
    R <- ncol(parameters_boot)
    NM <- length(models)
    test_models <- NULL

    if (NM>0){
        for (mm in 1L:NM){
            model_mm <- models[[mm]]
            pp <- names(models)[mm]
            ind_pp <- which( paste(parameters$par)==pp )
            y <- parameters[ind_pp,'est']
            dat <- data.frame(m=mod$moderator.grid, y=y, w=w)
            mod11 <- stats::lm(formula=model_mm, data=dat, weights=w)
            coef11 <- coef(mod11)
            parameters[ind_pp,'est'] <- predict(mod11)

            NC <- length(coef11)
            if (! bmod_missing){
                est_boot <- matrix(NA, nrow=NC, ncol=R)
                rr <- 1
                for (rr in 1:R){
                    dat$y <- parameters_boot[ind_pp,rr]
                    dat$w <- bmod$moderator_density_boot[,rr]
                    mod12 <- stats::lm(formula=model_mm, data=dat)
                    parameters_boot[ind_pp,rr] <- predict(mod12)
                    est_boot[,rr] <- coef(mod12)
                }
            }
            dfr1 <- data.frame(par=pp, coef=names(coef11), est=coef11)
            if (! bmod_missing){
                dfr1$se <- apply(est_boot, 1, stats::sd)
                dfr1$t <- dfr1$est / dfr1$se
                dfr1$p <- 2*stats::pnorm( -abs(dfr1$t) )
                # global Wald test for all parameters without intercept
                V <- stats::cov(t(est_boot))
                A <- matrix(0, nrow=NC-1, ncol=NC)
                for (cc in 1:(NC-1)){
                    A[cc,cc+1] <- 1
                }
                res <- lsem_wald_test(theta=coef11, V=V, A=A)
                dfr1[1,'chisq'] <- res$chisq
                dfr1[1,'df'] <- res$df
                dfr1[1,'p_wald'] <- res$p
            }
            rownames(dfr1) <- NULL
            test_models <- rbind(test_models, dfr1)
        }
    }

    #--- output
    res <- list(wald_test_global=wald_test_global, test_models=test_models,
                    parameters=parameters, parameters_boot=parameters_boot)
    return(res)
}
