## File Name: sirt_summary_print_call.R
## File Version: 0.01

sirt_summary_print_call <- function(CALL)
{
	cat("Call:\n", paste(deparse(CALL), sep = "\n", collapse = "\n"), 
				"\n\n", sep = "")
}				
