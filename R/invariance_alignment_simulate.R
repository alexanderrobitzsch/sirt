## File Name: invariance_alignment_simulate.R
## File Version: 0.108

invariance_alignment_simulate <- function(nu, lambda, err_var, mu, sigma, N,
    output="data", groupwise=FALSE)
{
    if (length(N)==1 & nrow(nu)>1){
        N <- rep(N, nrow(nu))
    }
    N_tot <- sum(N)
    G <- nrow(nu)
    I <- ncol(nu)
    items <- colnames(nu)
    if (is.null(items)){
        items <- paste0("I",1:I)
    }
    n_end <- cumsum(N)
    n_start <- c(1,n_end+1)[-c(G+1)]
    #* simulate data
    if (N[1]<Inf){
        group <- rep(1:G, N)
        dat <- matrix(NA, nrow=N_tot, ncol=I+1)
        colnames(dat) <- c("group",items)
        dat <- as.data.frame(dat)
        dat$group <- group
        for (gg in 1:G){
            N_gg <- N[gg]
            ind_gg <- seq(n_start[gg], n_end[gg])
            fac <- stats::rnorm(N_gg, mean=mu[gg], sd=sigma[gg])
            for (ii in 1:I){
                err_ii <- stats::rnorm(N_gg, mean=0, sd=sqrt(err_var[gg,ii]) )
                dat[ind_gg, ii+1] <- nu[gg,ii] + lambda[gg,ii]*fac + err_ii
            }
        }
        #-- output
        res <- dat
        if (output=="suffstat"){
            res <- list(mu=list(), Sigma=list(), N=list() )
            for (gg in 1:G){
                ind_gg <- which(dat$group==gg)
                res$mu[[gg]] <- colMeans(dat[ ind_gg, -1])
                res$Sigma[[gg]] <- stats::cov.wt(dat[ ind_gg, -1], method="ML")$cov
                res$N[[gg]] <- N[gg]
            }
        }
    }
    #*** only compute covariance matrices
    if (N[1]==Inf){
        res <- list(mu=list(), Sigma=list(), N=list() )
        for (gg in 1:G){
            lam_gg <- lambda[gg,]
            sig2_gg <- sigma[gg]^2
            res$mu[[gg]] <- nu[gg,] + lam_gg*mu[gg]
            lam_gg <- matrix(lam_gg, ncol=1)
            res$Sigma[[gg]] <- lam_gg %*% sig2_gg %*% t(lam_gg) + diag(err_var[gg,])
            res$N[[gg]] <- N[gg]
        }
        output <- "suffstat"
    }
    #*** group-wise output
    if ( (output=="suffstat") & groupwise ){
        res1 <- res
        res <- as.list(1:G)
        for (gg in 1:G){
            res[[gg]] <- list(mu=res1$mu[[gg]], Sigma=res1$Sigma[[gg]], N=res1$N[[gg]])
        }
    }

    #--- output
    return(res)
}
