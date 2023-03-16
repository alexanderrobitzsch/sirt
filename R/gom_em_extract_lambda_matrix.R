## File Name: gom_em_extract_lambda_matrix.R
## File Version: 0.01
## File Last Change: 2019-05-13

gom_em_extract_lambda_matrix <- function(lambda_logit, I, K)
{
    lambda <- matrix( stats::plogis(lambda_logit), nrow=I, ncol=K)
    return(lambda)
}
