## File Name: rasch_mml2_calcprob_missing1.R
## File Version: 1.092
## File Last Change: 2021-05-08



#**** calculation of probabilities in the missing data IRT model
rasch_mml2_calcprob_missing1 <- function( theta.k, b, beta, delta.miss, pjk,
        fixed.a=NULL, return_nonresponse_probs=FALSE, irtmodel="missing1")
{
    I <- length(b)
    if (is.null(fixed.a)){
        fixed.a <- rep(1,I)
    }
    TP <- nrow(theta.k)

    thetaM1 <- sirt_matrix2(theta.k[,1], nrow=I)
    thetaM2 <- sirt_matrix2(theta.k[,2], nrow=I)
    ab <- fixed.a*b
    M1 <- stats::plogis(fixed.a*thetaM1-ab)  # probability for correct item response

    if (irtmodel=="missing1"){
        # probability of a response for incorrect item responses
        M2a <- stats::plogis(thetaM2-beta)
        # probability of a response for correct item responses
        M2b <- stats::plogis(thetaM2-beta-delta.miss)
    } else {
        M2a <- stats::plogis(thetaM2-beta+delta.miss)
        M2b <- stats::plogis(thetaM2-beta)
    }
    pjk[,1,] <- 1 - M1
    pjk[,2,] <- M1
    # compute joint probability (R,Y)
    pjk[,1,] <- M2a*pjk[,1,]
    pjk[,2,] <- M2b*pjk[,2,]
    # P(R=0)=P(R=0|Y=0)P(Y=0)+P(R=0|Y=1)P(Y=1)
    m3 <- pjk[,3,] <- (1-M2a)*(1-M1) + (1-M2b)*M1
    if (return_nonresponse_probs){
        pjk[,1,] <- (1-M2a)*(1-M1) / m3
        pjk[,2,] <- (1-M2b)*M1 / m3
        pjk <- pjk[,1:2,]
    }
    return(pjk)
}


.calcprob.missing1 <- rasch_mml2_calcprob_missing1
