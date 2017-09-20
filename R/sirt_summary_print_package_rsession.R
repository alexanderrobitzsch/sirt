## File Name: sirt_summary_print_package_rsession.R
## File Version: 0.01
## File Last Change: 2017-09-20 09:50:10

sirt_summary_print_package_rsession <- function(pack)
{
    res <- TAM::tam_packageinfo(pack=pack)			
	cat(res,"\n")
	res <- TAM::tam_rsessinfo()
	cat(res, "\n\n")
}
