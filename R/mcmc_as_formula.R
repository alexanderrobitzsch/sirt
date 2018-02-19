## File Name: mcmc_as_formula.R
## File Version: 0.06

mcmc_as_formula <- function( string )
{
	string <- paste0( string , collapse = " " )
	string <- gsub("___ " , "___" , string, fixed=TRUE )
	form <- stats::as.formula(string)
	return(form)
}
