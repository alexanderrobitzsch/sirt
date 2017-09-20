## File Name: summary.xxirt.R
## File Version: 0.12
## File Last Change: 2017-09-20 10:47:26
#*******************************************************
# Summary for xxirt object
summary.xxirt <- function( object , digits = 3 , file=NULL , ...){

    # open sink
    sirt_osink( file = file )

	cat("-----------------------------------------------------------------\n")
    d1 <- utils::packageDescription("sirt")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
	cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
	cat("Computation Time:" , print(object$s2 - object$s1), "\n\n")
	
	cat("Call:\n", paste(deparse(object$CALL), sep = "\n", collapse = "\n"), 
				"\n\n", sep = "")	
					
#	modeltype <- object$irtmodel
		cat( "   " , object$ic$n , "Cases, " , object$ic$I , "Items, " , 
		        object$G , "Group(s)", # "," ,
				"\n")  

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
												
    cat( "AIC  = " , round( object$ic$AIC , 0 ) , " | penalty =" , round( object$ic$AIC - object$ic$deviance ,2 ) , 
			"   | AIC = -2*LL + 2*p  \n" )    
    cat( "AICc = " , round( object$ic$AICc , 0 ) ," | penalty =" , round( object$ic$AICc - object$ic$deviance ,2 ) )
		cat("    | AICc = -2*LL + 2*p + 2*p*(p+1)/(n-p-1)  (bias corrected AIC)\n" )   	
    cat( "BIC  = " , round( object$ic$BIC , 0 ) , " | penalty =" , round( object$ic$BIC - object$ic$deviance ,2 ) , 
			"   | BIC = -2*LL + log(n)*p  \n" )  
    cat( "CAIC = " , round( object$ic$CAIC , 0 ) ," | penalty =" , round( object$ic$CAIC - object$ic$deviance ,2 ) )
		cat("   | CAIC = -2*LL + [log(n)+1]*p  (consistent AIC)\n\n" )   

    cat("-----------------------------------------------------------------\n")
	cat("Trait Parameters\n")
	v1 <- object$customTheta$par
	print( round( v1 , digits )  )
    cat("-----------------------------------------------------------------\n")
	cat("Item Parameters \n")
	obji <- object$par_items_summary
	obji[,-c(1:2)] <- round( obji[,-c(1:2)] , digits)
	print( obji ) 

	# close sink
    sirt_csink( file = file )		
}
#*******************************************************
