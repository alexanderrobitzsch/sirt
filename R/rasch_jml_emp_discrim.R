## File Name: rasch_jml_emp_discrim.R
## File Version: 0.08
## File Last Change: 2018-12-30



# Function for calculating empirical discrimination
# slope estimation (WINSTEPS manual p. 300)
rasch_jml_emp_discrim <- function( theta, b, dat, dat.resp=1-is.na(dat.resp), freq )
{
    N <- length(theta)
    I <- length(b)
    pni <- .prob.rasch( theta=theta, b=b )
    bM <- matrix(b, nrow=N, ncol=I, byrow=TRUE)
    thetaM <- matrix(theta, nrow=N, ncol=I)
    tbdiff <- thetaM - bM
    tdf <- tbdiff * dat.resp * freq
    t1 <- colSums( ( dat - pni ) * tdf )
    t2 <- colSums( pni * ( 1 - pni ) * tbdiff * tdf)
    res <- 1 + t1/t2
    return(res)
}

# emp.discr <- rasch_jml_emp_discrim
