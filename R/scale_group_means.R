## File Name: scale_group_means.R
## File Version: 0.08
## File Last Change: 2021-06-02

scale_group_means <- function(M, SD, probs=NULL, M_target=0, SD_target=1)
{
    if (is.null(probs)){
        probs <- rep(1,length(M))
    }
    probs <- probs/sum(probs)
    M1 <- sum(probs*M)
    SD1 <- sqrt( sum(probs*SD^2)+sum(probs*(M-M1)^2) )
    M_z <- ( M - M1 ) / SD1
    SD_z <- SD / SD1
    M1_trafo <- M_target + SD_target*M_z
    SD1_trafo <- SD_target*SD_z
    model_M <- stats::lm(M1_trafo~M)
    # trafo y=a*x+b
    a <- coef(model_M)[2]
    b <- coef(model_M)[1]
    #--- output
    res <- list(M1=M1, SD1=SD1, M_z=M_z, SD_z=SD_z, M_trafo=M1_trafo,
                    SD_trafo=SD1_trafo, a=a, b=b )
    return(res)
}
