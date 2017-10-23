## File Name: rm_facets_print_progress_deviance.R
## File Version: 0.01
## File Last Change: 2017-10-02 15:45:14


rm_facets_print_progress_deviance <- function( dev, dev0, digits_deviance, iter )
{
	cat( paste( "   Deviance = "  , round( dev , digits_deviance ) , 
		if (iter > 1 ){ " | Deviance change = " } else {""} ,
		if( iter>1){ round( - dev + dev0 , digits_deviance + 2)} else { ""}	,"\n",sep="") )
}
