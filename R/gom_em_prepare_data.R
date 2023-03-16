## File Name: gom_em_prepare_data.R
## File Version: 0.05
## File Last Change: 2019-05-18

gom_em_prepare_data <- function(dat, weights, model)
{
    dat0 <- dat
    dat.resp <- 1-is.na(dat)
    dat[ is.na(dat) ] <- 0
    N <- nrow(dat)
    I <- ncol(dat)
    if (is.null(weights)){
        weights <- rep(1,N)
    } else {
        if (model=="GOMRasch"){
            stop("'GOMRasch' cannot handle weights!\n")
        }
    }
    dat2 <- as.matrix(dat)
    dat2.resp <- as.matrix(dat.resp)
    # indicator matrix
    dat2.ind0 <- dat2.resp*(dat2==0)
    dat2.ind1 <- dat2.resp*(dat2==1)
    dat2.ind <- as.matrix( cbind( dat2.ind0, dat2.ind1 ) )

    #--- output
    res <- list( dat0=dat0, dat2=dat2, dat2.resp=dat2.resp, dat2.ind=dat2.ind, N=N,
                    I=I, weights=weights )
    return(res)
}
