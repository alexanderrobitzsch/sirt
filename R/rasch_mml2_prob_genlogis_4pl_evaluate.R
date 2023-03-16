## File Name: rasch_mml2_prob_genlogis_4pl_evaluate.R
## File Version: 0.08
## File Last Change: 2019-05-10


rasch_mml2_prob_genlogis_4pl_evaluate <- function(theta.k, b, fixed.a, fixed.c,
    fixed.d, alpha1, alpha2, Qmatrix, h, h1, h2, incr=NULL)
{
    if (!is.null(incr) & (h!=0) ){
        for (vv in incr){
            if (vv=="b"){ b <- b + h }
            if (vv=="a"){ fixed.a <- fixed.a + h }
            if (vv=="c"){ fixed.c <- fixed.c + h }
            if (vv=="d"){ fixed.d <- fixed.d + h }
            if (vv=="alpha1"){ alpha1 <- alpha1 + h }
            if (vv=="alpha2"){ alpha2 <- alpha2 + h }
        }
    }
    pjk <- prob_genlogis_4pl(theta=theta.k, b=b, a=fixed.a, c=fixed.c,
                d=fixed.d, alpha1=alpha1, alpha2=alpha2, Qmatrix=Qmatrix)
    pjk <- ( pjk + h1 ) / h2
    pjk.M <- t(pjk)
    qjk.M <- 1 - pjk.M
    #-- output
    res <- list(pjk.M=pjk.M, qjk.M=qjk.M)
    return(res)
}
