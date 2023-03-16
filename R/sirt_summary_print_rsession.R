## File Name: sirt_summary_print_rsession.R
## File Version: 0.04


sirt_summary_print_rsession <- function()
{
    res <- TAM::tam_rsessinfo()
    cat(res, '\n\n')
}
