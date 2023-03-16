## File Name: mcmc_rename_define_symbols.R
## File Version: 0.05
## File Last Change: 2018-12-30

mcmc_rename_define_symbols <- function()
{
    trans <- c("X", "Z", "M")
    orig <- c("[", "]", ",")
    res <- list(trans=trans, orig=orig)
    return(res)
}
