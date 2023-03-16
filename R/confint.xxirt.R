## File Name: confint.xxirt.R
## File Version: 0.11
## File Last Change: 2018-12-30

###########################################
# confidence interval
confint.xxirt <- function( object, parm, level=.95,  ... )
{
    c1 <- coef.xxirt(object)
    v1 <- vcov.xxirt(object)
    if ( ! missing(parm) ){
        c1 <- c1[parm]
        v1 <- v1[ parm, parm]
    }

    q1 <- ( 1 - level ) / 2
    q2 <- 1 - ( 1 - level ) / 2
    quant <- stats::qnorm(q2)
    se <- sqrt( diag(v1) )
    res <- data.frame( "a"=c1 - quant*se, "b"=c1 + quant*se )
    colnames(res)[1] <- paste0( round( 100*q1,1 ), " %")
    colnames(res)[2] <- paste0( round( 100*q2,1 ), " %")
    rownames(res) <- names(c1)
    return(res)
}
###############################################
