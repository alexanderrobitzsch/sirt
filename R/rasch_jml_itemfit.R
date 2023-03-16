## File Name: rasch_jml_itemfit.R
## File Version: 0.05



#*** item fit Rasch model
rasch_jml_itemfit <- function( theta0, b, dat )
{
    dat9 <- dat
    dat9[ is.na(dat)] <- 9
    ind <- is.finite(theta0)
    dat2 <- dat9[ ind, ]
    dat2.resp <- 1 * ( dat2 !=9 )
    theta0 <- theta0[ ind ]
    p <- .prob.rasch( theta0, b )
    v <- p * (1 - p)
    z2 <- dat2.resp * (dat2 - p)/v
    z2d <- z2 * dat2.resp
    infit <- colSums( z2d * v ) / colSums( v * dat2.resp )
    outfit <- colSums(z2d) / colSums( dat2.resp )
    res <- data.frame(infit, outfit )
    return(res)
}

rasch.itemfit <- rasch_jml_itemfit

