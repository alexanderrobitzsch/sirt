## File Name: xxirt_compute_posterior.R
## File Version: 0.301


##-- xxirt: compute posterior
xxirt_compute_posterior <- function( prior_Theta, p.xi.aj, group,
                G, weights, dat1, dat_resp, maxK, group_index, dat1_resp,
                eps=1e-100, customTheta=NULL)
{
    N <- nrow(dat_resp)
    TP <- ncol(p.xi.aj)
    I <- ncol(dat1)
    # posterior distribution
    if (customTheta$person_covariates){
        prior1 <- t(prior_Theta)
    } else {
        prior1 <- t( prior_Theta[, group ] )
    }

    # p.xi.aj[ is.na(p.xi.aj) ] <- eps
    p.aj.xi <- prior1 * p.xi.aj
    p1 <- p.aj.xi
    p.aj.xi <- p.aj.xi / rowSums( p.aj.xi )

    # expected counts
    n.ik <- array( 0, dim=c(I,maxK, TP,G) )
    N.ik <- array( 0, dim=c(I,maxK, TP) )
    pi.k <- matrix( 0, nrow=TP, ncol=G )
    for (gg in 1L:G){
        ind_gg <- group_index[[gg]]
        p.aj.xi.gg <- as.matrix( p.aj.xi[ind_gg, ] )
        dat1_resp_gg <- dat1_resp[ ind_gg,, ]
        weights_gg <- weights[ind_gg]
        for (kk in 1L:maxK){
            dat1_resp_gg_kk <- as.matrix(dat1_resp_gg[,,kk])
            n.ik[,kk,,gg] <- sirt_rcpp_xxirt_compute_posterior_expected_counts(
                                    dat1_resp_gg=dat1_resp_gg_kk, p_aj_xi_gg=p.aj.xi.gg,
                                    weights_gg=weights_gg)
        }
        N.ik <- N.ik + n.ik[,,,gg]
        pi.k[,gg] <- colSums( p.aj.xi.gg * weights[ ind_gg ] )
    }  # end gg
    res <- list( p.aj.xi=p.aj.xi, n.ik=n.ik, N.ik=N.ik, N.k=pi.k, post_unnorm=p1 )
    return(res)
}
