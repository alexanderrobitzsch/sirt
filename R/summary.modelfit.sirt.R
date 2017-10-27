## File Name: summary.modelfit.sirt.R
## File Version: 1.03


##############################################################
# summary modelfit.sirt
summary.modelfit.sirt <- function( object , ... )
{	
	cat("Test of Global Model Fit\n")
	sirt_summary_print_objects(obji=object$modelfit.test, digits=5, from=2)
	
	cat("\nFit Statistics\n")
	sirt_summary_print_objects(obji=object$modelfit.stat, digits=5, from=1,
				rownames_null=FALSE)
}
#################################################################	


