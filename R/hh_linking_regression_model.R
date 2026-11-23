## File Name: hh_linking_regression_model.R
## File Version: 0.158

hh_linking_regression_model <- function(modlist, partype, verbose=TRUE,
        prec_weights="none", rob_weights="none", c=NULL, iter_max=500, tol=1e-5,
        eps=0.001)
{
    parinfo <- modlist$parinfo
    X <- modlist$X
    G <- modlist$G
    I <- modlist$I

    ind <- which( parinfo$partype==partype )

    y <- parinfo[ ind, 'est']

    V <- modlist$varmat
    V <- V[ind, ind]

    if (partype=='a'){
        A <- diag(1/ y)
        V <- A %*% V %*% t(A)
    }

    if (partype=='a'){
        y <- pmax(y, 0.05)
        y <- log( abs(y) )
    }
    if (partype=='d'){
        sigma <- c(1, modlist$sigma)
        parinfo1 <- parinfo[ parinfo$partype=='a', ]
        for (gg in 2L:G){
            ind_gg <- which(parinfo1$group==gg)
            aig <- parinfo1[ind_gg, 'est']
            X[ind_gg,gg-1] <- aig / sigma[gg]
        }
    }

    #-- include pairwise modeling
    P <- modlist$P
    Z <- modlist$Z
    X <- ( P %*% X )[, 1L:(G-1) ]
    y <- as.vector( P %*% y )
    V <- P %*% V %*% t(P)
    NP <- length(y)

    weights <- hh_linking_apply_precision_weights(V=V, prec_weights=prec_weights)

    if (verbose){
        if (partype=='a'){
            cat('-- Estimate sigma\n')
        } else {
            cat('-- Estimate mu\n')
        }
        utils::flush.console()
    }

    coef0 <- rep(0, G-1)
    iter <- 1
    iterate <- TRUE

    iter_max <- 50

    #- fit regression model
    while (iterate){

        mod <- stats::lm( y ~ 0 + X, weights=weights)
        coef1 <- mod$coef
        change <- max(abs(coef0-coef1))
        if (verbose){
            cat( paste0(' * Iteration ', iter, ': Parameter change=',
                        round(change,8), '\n' ))
            utils::flush.console()
        }
        coef0 <- coef1
        fitted_values <- mod$fitted.values
        residuals <- y - fitted_values

        #- apply weighting
        weights <- hh_linking_apply_precision_weights(V=V, prec_weights=prec_weights)
        robw <- hh_linking_apply_robust_weights( residuals=residuals,
                        rob_weights=rob_weights, V=V, eps=eps, c=c )
        weights <- weights*robw
        weights <- normalize_weights(x=weights)
        W <- diag(weights)
        iter <- iter + 1
        if ((change < tol) | (iter > iter_max) ) {
            iterate <- FALSE
        }
    }

    #- output
    delta <- mod$coef[1L:(G-1)]
    names(delta) <- paste0('G',2L:G)

    if (partype=='a'){
        fitted_values <- exp(fitted_values)
        delta <- exp(delta)
        modlist$sigma <- delta
    }

    #- output
    res <- list(delta=delta, partype=partype, parinfo=parinfo, iter=iter,
                weights=weights, robw=robw, y=y, W=W, X=X, prec_weights=prec_weights )
    return(res)
}
