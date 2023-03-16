## File Name: predict.btm.R
## File Version: 0.07
## File Last Change: 2020-04-18

predict.btm <- function(object, data=NULL, ...)
{
    dat <- data
    if (is.null(dat)){
        dat <- object$dat
    }
    effects <- object$effects
    teams <- paste(effects$individual)
    pars <- object$pars
    ignore.ties <- object$ignore.ties
    delta <- pars[ pars$par=="delta", "est"]
    eta <- pars[ pars$par=="eta", "est"]
    wgt.ties <- object$wgt.ties
    #- identify teams
    ind1 <- match( paste(dat[,1]), teams)
    ind2 <- match( paste(dat[,2]), teams)

    #- compute probabilities
    theta <- effects$theta
    ND <- nrow(dat)
    M1 <- matrix(0, nrow=ND, ncol=3)
    M1[,1] <- theta[ ind1 ] + eta
    M1[,2] <- theta[ ind2 ]
    M1[,3] <- delta + ( theta[ ind1] + theta[ ind2 ] + eta ) * wgt.ties
    probs <- exp(M1)
    probs <- probs / rowSums(probs)
    #-- output
    return(probs)
}
