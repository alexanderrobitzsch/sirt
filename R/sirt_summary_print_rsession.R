## File Name: sirt_summary_print_rsession.R
## File Version: 0.01
## File Last Change: 2017-10-02 11:15:49


sirt_summary_print_rsession <- function()
{
	res <- TAM::tam_rsessinfo()
	cat(res, "\n\n")
}
