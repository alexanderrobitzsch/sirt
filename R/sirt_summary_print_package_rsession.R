## File Name: sirt_summary_print_package_rsession.R
## File Version: 0.04

sirt_summary_print_package_rsession <- function(pack)
{
    res <- TAM::tam_packageinfo(pack=pack)
    cat(res,"\n")
    res <- TAM::tam_rsessinfo()
    cat(res, "\n")
}
