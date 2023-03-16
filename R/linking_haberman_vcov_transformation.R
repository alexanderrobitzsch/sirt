## File Name: linking_haberman_vcov_transformation.R
## File Version: 0.12
## File Last Change: 2021-12-09

linking_haberman_vcov_transformation <- function( H1, aj_vcov )
{
    if (is.null(aj_vcov)){
        aj_se <- NA
    } else {
        aj_vcov <- H1 %*% aj_vcov %*% t(H1)
        aj_se <- c( sqrt( diag( aj_vcov ) ) )
    }
    res <- list( vcov=aj_vcov, se=aj_se )
    return(res)
}
