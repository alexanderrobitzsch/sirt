## File Name: summary.lsem.permutationTest.R
## File Version: 0.13
## File Last Change: 2017-01-31 18:33:40
######################################################
summary.lsem.permutationTest <- function( object , file=NULL , digits=3 , ... ){

	# open sink for a file
	CDM::osink( file=file , suffix="__SUMMARY.Rout" )

	cat("-----------------------------------------------------------------\n")
	cat("Permutation Test for Local Structural Equation Model \n\n")
		
	cat( package_version_date("sirt") , "\n" )
	cat( package_version_date("lavaan") , "\n" )	
	cat( package_version_date("lavaan.survey") , "\n" )					
	
	cat("\nFunction 'sirt::lsem.permutationTest' \n\n")
	
	cat("Call:\n", paste(deparse(object$CALL), sep = "\n", collapse = "\n"), 
				"\n\n", sep = "")	
	
	cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
	cat("Computation Time:" , print(object$s2 - object$s1), "\n\n")
    
	cat( "Number of permutations = " , object$B , "\n")
	cat( paste0( "Number of observations = " , round(object$N,digits) ) , "\n")
	cat( paste0( "Bandwidth factor = " , round(object$h,digits) ) , "\n")
	cat( paste0( "Bandwidth = " , round(object$bw,digits) ) , "\n")
    cat( paste0( "Number of focal points for moderator = " , 
						length(object$moderator.grid ) ) , "\n")	
	
	cat("\nlavaan Model\n")
	cat(object$lavmodel)
	
	cat("\n\n")
	
	cat("Global Test Statistics\n\n")
	obji <- object$teststat		
	VV <- ncol(obji)
	for (vv in 2:VV){
		obji[,vv] <- round( obji[,vv] , digits )
					}					
	print(obji)			
	cat("\n")
	cat("Pointwise Test Statistics\n\n")
	obji <- object$parameters_pointwise_test		
	vars <- c("est" , "p")
	for (vv in vars){
		obji[,vv] <- round( obji[,vv] , digits )
					}
	obji <- obji[ , c("par" , "parindex" , "moderator" , "est" , "p") ]					
	rownames(obji) <- NULL					
	print(obji)	

	# close file
	CDM::csink(file)
	
}
######################################################			
