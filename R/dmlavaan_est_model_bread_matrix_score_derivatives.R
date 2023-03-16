## File Name: dmlavaan_est_model_bread_matrix_score_derivatives.R
## File Version: 0.112


dmlavaan_est_model_bread_matrix_score_derivatives <- function(fun, args,
        partable, scores, symmetric=TRUE, h=1e-4, remove_duplicated=TRUE,
        parnames=NULL)
{
    N <- nrow(scores)
    NPU <- max(partable$pnid)
    args$do.fit <- FALSE
    hess <- array(NA, dim=c(N, NPU, NPU) )
    partable1 <- partable
    partable1$start <- partable$est
    for (pp in 1:NPU){
        partable1a <- partable1
        ind_pp <- which( partable$pnid==pp )
        val <- partable1$est[ind_pp[1]]
        h1 <- ifelse( abs(val) > 1, abs(val)*h, h )
        #* add increment
        partable1a$start[ ind_pp ] <- partable1$start[ ind_pp ] + h1
        partable1a$est <- partable1a$start
        args$model <- partable1a
        mod2 <- do.call(what=fun, args=args)
        scores2 <- lavaan::lavScores(object=mod2, remove.duplicated=remove_duplicated)
        #* substract increment
        if (symmetric){
            partable1a$start[ ind_pp ] <- partable1$start[ ind_pp ] - h1
            partable1a$est <- partable1a$start
            args$model <- partable1a
            mod2 <- do.call(what=fun, args=args)
            scores3 <- lavaan::lavScores(object=mod2, remove.duplicated=remove_duplicated)
            fac <- 2
        } else {
            scores3 <- scores
            fac <- 1
        }
        hess[,pp,] <- ( scores2-scores3 ) / (fac*h1)
    }
    B <- matrix(0, nrow=NPU, ncol=NPU)
    pp <- 1; hh <- 1
    for (pp in 1:NPU){
        for (hh in 1:NPU){
            B[pp,hh] <- mean( hess[,pp,hh] )
        }
    }
    B <- -B
    # B <- ( B + t(B) ) / 2
    if (!is.null(parnames)){
        colnames(B) <- rownames(B) <- parnames
    }
    #--- output
    return(B)
}
