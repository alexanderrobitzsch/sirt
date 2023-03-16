## File Name: sirt_pmvnorm.R
## File Version: 0.02
## File Last Change: 2023-03-08


sirt_pmvnorm <- function(... )
{
    TAM::require_namespace_msg('mvtnorm')
    res <- mvtnorm::pmvnorm(...)
    return(res)
}
