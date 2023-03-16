## File Name: predict_scale_group_means.R
## File Version: 0.03
## File Last Change: 2021-09-22


predict_scale_group_means <- function(object, M, SD)
{
    a <- object$a
    b <- object$b
    M_trafo <- a*M+b
    SD_trafo <- a*SD
    #-- output
    res <- list(M_trafo=M_trafo, SD_trafo=SD_trafo)
    return(res)
}
