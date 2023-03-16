## File Name: rasch_mml2_impute_data.R
## File Version: 0.127
## File Last Change: 2021-05-08


rasch_mml2_impute_data <- function(e1, b, beta, delta.miss, pjk, fixed.a,
                            dat, dat.resp, irtmodel, nimps=10)
{
    theta.k <- e1$theta.k
    if (is.vector(theta.k)){
        theta.k <- matrix(theta.k, ncol=1)
    }

    TP <- nrow(theta.k)
    N <- nrow(dat.resp)
    I <- ncol(dat.resp)
    #--- missing1
    if (irtmodel=="missing1"){
        # P(YR|theta,xsi) is given
        # P(Y|R=0,theta, xsi)=P(Y=1,R=0,theta,xsi) /
        #                         { P(Y=0,R=0, theta,xsi)+P(Y=1,R=0, theta,xsi) }
        probs <- rasch_mml2_calcprob_missing1( theta.k=theta.k, b=b, beta=beta,
                        delta.miss=delta.miss, pjk=pjk, fixed.a=fixed.a,
                        return_nonresponse_probs=TRUE, irtmodel=irtmodel)
    }

    #--- raschtype
    if (irtmodel=="raschtype"){
        probs <- array(NA, dim=c(I,2,TP))
        probs[,1,] <- 1-t(pjk)
        probs[,2,] <- t(pjk)
    }
    post <- rowCumsums.sirt(matr=e1$f.qk.yi)
    dat_implist <- list()
    for (uu in 1:nimps){
        dat0 <- dat
        rn <- stats::runif(N)
        indices <- rowIntervalIndex.sirt(matr=post,rn=rn)
        for (ii in 1:I){
            probs_ii <- t(probs[ii,,indices])
            rn_ii <- stats::runif(N)
            y_lat <- 1*(rn_ii>probs_ii[,1])
            ind <- which(paste(dat[,ii])==2)
            if (length(ind)>0){
                dat0[ind,ii] <- y_lat[ind]
            }
        }  # end ii (items)
        dat_implist[[uu]] <- dat0
    }  # end uu (number of imputations)
    #--- output
    return(dat_implist)
}
