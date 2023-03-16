## File Name: sirt_var.R
## File Version: 0.01
## File Last Change: 2022-04-10

sirt_var <- function(x, method="ML", na.rm=TRUE)
{
    v1 <- stats::var(x, na.rm=na.rm)
    N <- sum(1-is.na(x))
    v1 <- v1 * (N-1) / N
    return(v1)

}
