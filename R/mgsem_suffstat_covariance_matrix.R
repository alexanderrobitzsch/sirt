## File Name: mgsem_suffstat_covariance_matrix.R
## File Version: 0.099

mgsem_suffstat_covariance_matrix <- function(suffstat)
{
    G <- length(suffstat)
    dfr <- NULL
    # gg <- 1

    for (gg in 1:G){
        suffstat_gg <- suffstat[[gg]]
        names_gg <- names(suffstat_gg)
        mu <- suffstat_gg[[ intersect( c('mu','M'), names_gg) ]]
        Sigma <- suffstat_gg[[ intersect( c('S','Sigma'), names_gg) ]]
        N <- suffstat_gg$N
        vars <- names(mu)
        I1 <- length(mu)
        if (is.null(vars)){
            vars <- paste0('V',1:I1)
        }
        #** mu
        label <- paste0('mu_G', gg, '_', vars)
        se <- sqrt( diag(Sigma) / N )
        dfr1 <- data.frame(type='mu', group=gg, index1=1:I1, index2=NA, label=label,
                    par=mu, se=se)
        V1 <- Sigma / N
        rownames(V1) <- colnames(V1) <- label

        #** sigma
        c1 <- rbind( 1:I1, 1:I1)
        c2 <- utils::combn(x=I1, m=2)
        c3 <- rbind( t(c1), t(c2))
        c3 <- c3[ order(c3[,1]), ]
        label <- paste0('Sigma_G', gg, '_', vars[c3[,1]], '_', vars[c3[,2]] )
        par <- Sigma[ c3[,1:2] ]

        #* compute duplication matrix
        K <- mgsem_duplication_matrix(x=Sigma)

        dfr2 <- data.frame(type='sigma', group=gg, index1=c3[,1], index2=c3[,2],
                        label=label, par=par)
        N2 <- nrow(dfr2)
        V2 <- matrix(0, nrow=N2, ncol=N2)
        V2 <- 2* ( K %*% ( kronecker(X=Sigma, Y=Sigma) ) %*% t(K) ) / N
        rownames(V2) <- colnames(V2) <- label
        dfr2$se <- sqrt(diag(V2))
        V3 <- mgsem_bdiag(x1=V1, x2=V2)
        dfr <- rbind(dfr, dfr1, dfr2)
        if (gg==1){
            V <- V3
        } else {
            V <- mgsem_bdiag(x1=V, x2=V3)
        }

    }    # end gg

    #--- output
    res <- list(suffstat_pars=dfr, suffstat_vcov=V)
    return(res)
}
