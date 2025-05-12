## File Name: xxirt_mstep_ThetaParameters.R
## File Version: 0.249


xxirt_mstep_ThetaParameters <- function( customTheta, G, eps,
            mstep_iter, N.k, par1, mstep_reltol, Theta, penalty_fun_theta=NULL,
            weights=NULL, p.aj.xi=NULL)
{
    like_Theta <- function( x, ... )
    {
        par1 <- customTheta$par
        par1[ customTheta$est ] <- x
        arg_list <- list( par=par1, Theta=Theta, G=G )
        if (customTheta$person_covariates){
            arg_list$X <- customTheta$X
            prior1 <- t( do.call( what=customTheta$P, args=arg_list )     )
            # multiply by individual log-likelihood
            ll2 <-  - sum( weights * p.aj.xi * log( prior1 + eps ) )
        } else {
            mod1 <- do.call( what=customTheta$P, args=arg_list )
            ll2 <- - sum( N.k * log( mod1 + eps ) )
        }

        NP <- length(customTheta$prior)
        pen <- 0
        if ( NP > 0 ){
            for (pp in 1L:NP){
                if ( ! is.na( customTheta$prior[pp] ) ) {
                    prior_pp <- customTheta$prior[pp]
                    # val <- par1[ names(customTheta$prior) ]
                    val <- par1[pp]
                    arg_list <- list( val,customTheta$prior_par1[pp],
                                                customTheta$prior_par2[pp]  )
                    prior_val <- do.call( prior_pp, arg_list )
                    pen <- pen - log( prior_val + eps )
                }  # end if prior
            }  # end pp
        }
        ll2 <- ll2 + pen
        if (!is.null(penalty_fun_theta)){
            pen <- penalty_fun_theta(x=x)
            ll2 <- ll2+pen
        }
        return(2*ll2)
    }
    #----- end definition likelihood function

    mstep_method <- 'BFGS'
    lower <- upper <- NULL
    if (customTheta$some_bound){
        mstep_method <- 'L-BFGS-B'
        lower <- customTheta$lower
        upper <- customTheta$upper
    }
    arg_control <- list(maxit=mstep_iter)
    if ( mstep_method %in% c( 'BFGS')){
        arg_control$reltol <- mstep_reltol
    }
    if ( mstep_method %in% c( 'L-BFGS-B')){
        arg_control$factr <- mstep_reltol
    }
    # method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
    arg_list <- list( par=par1, fn=like_Theta, method=mstep_method, control=arg_control)
    if (customTheta$some_bound){
        arg_list$lower <- lower
        arg_list$upper <- upper
    }
    mod <- do.call( what=stats::optim, args=arg_list )
    par1 <- mod$par
    ll2 <- mod$value
    customTheta$par[ customTheta$est ] <- par1
    par1 <- xxirt_ThetaDistribution_extract_freeParameters( customTheta=customTheta )

    pen_theta <- 0
    if (!is.null(penalty_fun_theta)){
        pen_theta <- penalty_fun_theta(x=par1)
    }

    #--- output
    res <- list( par1=par1, ll2=ll2, customTheta=customTheta, pen_theta=pen_theta)
    return(res)
}
