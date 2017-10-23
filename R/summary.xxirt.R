## File Name: summary.xxirt.R
## File Version: 0.16
## File Last Change: 2017-10-02 22:55:23
#*******************************************************
# Summary for xxirt object
summary.xxirt <- function( object , digits = 3 , file=NULL , ...)
{

    # open sink
    sirt_osink( file = file )

	cat("-----------------------------------------------------------------\n")
	#- package and R session
    sirt_summary_print_package_rsession(pack="sirt")	
	
	#- print call
	sirt_summary_print_call(CALL=object$CALL)
	
	#-- print computation time
	sirt_summary_print_computation_time_s1(object=object)
					
#	modeltype <- object$irtmodel
	cat( "   " , object$ic$n , "Cases, " , object$ic$I , "Items, " , 
		        object$G , "Group(s)", 	"\n")  

    cat("-----------------------------------------------------------------\n")
	cat( "Number of iterations =" , object$iter , "\n" )
    cat( "Deviance = " , round( object$deviance , 2 ) , " | " )
    cat( "Log Likelihood = " , round( -object$deviance/2 , 2 ) , "\n" )	
    cat( "Number of persons = " , object$ic$n , "\n" )    

    cat( "Number of estimated parameters = " , object$ic$np , "\n" )    
    cat( "  Number of estimated item parameters = " , object$ic$np.item , 
				"\n" )    	
    cat( "  Number of estimated distribution parameters = " , object$ic$np.Theta , 
				"\n\n" )    
												
	#--- information criteria
	rm_summary_information_criteria(object=object)

    cat("-----------------------------------------------------------------\n")
	cat("Trait Parameters\n")
	obji <- object$customTheta$par
	sirt_summary_print_objects(obji=obji, digits=digits, from=1)	
	
    cat("-----------------------------------------------------------------\n")
	cat("Item Parameters \n")
	obji <- object$par_items_summary
	sirt_summary_print_objects(obji=obji, digits=digits, from=2)	

	# close sink
    sirt_csink( file = file )		
}
#*******************************************************
