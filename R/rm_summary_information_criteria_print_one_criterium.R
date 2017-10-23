## File Name: rm_summary_information_criteria_print_one_criterium.R
## File Version: 0.03
## File Last Change: 2017-10-03 10:24:02

rm_summary_information_criteria_print_one_criterium <- function(ic, crit, desc, digits_crit=0, digits_penalty=2)
{
	val <- ic[[crit]]
	deviance <- ic[["deviance"]]
	crit_label0 <- crit
	if ( nchar(crit_label0) == 3){
		crit_label0 <- paste0( crit_label0 , " " )
	}	
	crit_label <- paste0( crit_label0 , " = " )
	penalty <- val - deviance
    cat( crit_label , round( val , digits_crit ) , " | penalty =" , round( penalty ,digits_penalty ) , 
			"   |" , desc , "  \n" ) 
}
