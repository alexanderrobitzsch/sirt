## File Name: rasch_mml2_estep_missing1.R
## File Version: 1.111



#*** E step missing data IRT model
rasch_mml2_estep_missing1 <- function( dat2, dat2.resp, theta.k, b, beta, delta.miss, I, CC,
        TP, group_, pi.k, pjk, weights, fixed.a=fixed.a, est.a=NULL, irtmodel="missing1" )
{
    eps <- 1e-16
    if (is.null(est.a)){
        est.a <- rep(0,I)
    }
a0 <- Sys.time()
    # probability correct response
    pjk <- rasch_mml2_calcprob_missing1( theta.k=theta.k, b=b, beta=beta,
                delta.miss=delta.miss, pjk=pjk, fixed.a=fixed.a, irtmodel=irtmodel )

    #** calculate likelihood
    probs_ <- as.matrix( array( pjk, dim=c(I,CC*TP) ) )
    f.yi.qk <- sirt_rcpp_probs_pcm_groups_C( dat=dat2, dat_resp=dat2.resp, group=group_,
                    probs=probs_, CC=CC, TP=TP )

    #*** calculate expected counts
    e1 <- sirt_rcpp_calccounts_pcm_groups_C( dat=dat2, dat_resp=dat2.resp, group=group_,
                fyiqk=f.yi.qk, pik=pi.k, CC=CC, weights=weights )
    e1$f.yi.qk <- f.yi.qk
    v1 <- array( e1$nik, dim=c(I,CC,TP) )
    e1$pjk <- pjk
    e1$n.k <- e1$count_pik
    e1$r.jk <- e1$n.jk <- NULL
    e1$n.ik <- v1
    e1$f.qk.yi <- e1$fqkyi
    e1$theta.k <- theta.k
    e1$irtmodel <- irtmodel
    return(e1)
}


# cat( " * calc counts") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1

.e.step.missing1 <- rasch_mml2_estep_missing1
