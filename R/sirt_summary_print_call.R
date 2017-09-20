## File Name: sirt_summary_print_call.R
## File Version: 0.01
## File Last Change: 2017-09-20 10:53:09

sirt_summary_print_call <- function(CALL)
{
	cat("Call:\n", paste(deparse(CALL), sep = "\n", collapse = "\n"), 
				"\n\n", sep = "")
}				
