## File Name: mcmc_rename_define_symbols.R
## File Version: 0.061

mcmc_rename_define_symbols <- function()
{
    trans <- c('X', 'Z', 'M')
    orig <- c('[', ']', ',')
    res <- list(trans=trans, orig=orig)
    return(res)
}
