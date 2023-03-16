## File Name: ruvn.R
## File Version: 0.02
## File Last Change: 2022-04-17

ruvn <- function(N, mean=0, sd=1, exact=TRUE)
{
    x <- stats::rnorm(N, mean=mean, sd=sd)
    y <- x
    if (exact){
        v1 <- sirt_var(x)
        m1 <- mean(x)
        y <- mean+(x-m1)/sqrt(v1)*sd
    }

    #--- output
    return(y)
}
