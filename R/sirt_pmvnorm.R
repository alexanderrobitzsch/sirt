## File Name: sirt_pmvnorm.R
## File Version: 0.01


sirt_pmvnorm <- function(... )
{
    TAM::require_namespace_msg("mvtnorm")
    res <- mvtnorm::pmvnorm(...)
    return(res)
}
