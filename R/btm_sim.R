## File Name: btm_sim.R
## File Version: 0.03
## File Last Change: 2019-03-10

btm_sim <- function(theta, eta=0, delta=-99, repeated=FALSE)
{
    N <- length(theta)
    dat <- t( utils::combn(N,2) )
    colnames(dat) <- c("id1", "id2")
    dat <- as.data.frame(dat)
    if (repeated){
        dat2 <- dat[, c("id2","id1") ]
        colnames(dat2) <- colnames(dat)
        dat <- rbind(dat, dat2 )
    }
    N1 <- nrow(dat)
    probs <- matrix(0, nrow=N1, ncol=3)
    probs[,1] <- eta + theta[ dat$id1 ]
    probs[,2] <- theta[ dat$id2 ]
    probs[,3] <- delta + ( eta + theta[ dat$id1 ] + theta[dat$id2 ] ) / 2
    probs <- exp(probs)
    probs <- probs / rowSums(probs)
    probs1 <- rowCumsums.sirt(matr=probs)
    rn <- stats::runif(N1)
    result <- rowIntervalIndex.sirt(matr=probs1,rn=rn)
    result <- 1*(result==1) + .5*(result==3) + 0*(result==2)
    dat <- data.frame( dat, result=result)
    rownames(dat) <- NULL
    return(dat)
}


