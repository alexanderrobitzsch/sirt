## File Name: summary.lsem.R
## File Version: 0.16
## File Last Change: 2017-01-31 18:34:11

#############################################
# summary lsem
summary.lsem <- function( object , file=NULL , digits=3 , ... ){

	# open sink for a file
	CDM::osink( file=file , suffix="__SUMMARY.Rout" )

	cat("-----------------------------------------------------------------\n")
	cat("Local Structural Equation Model \n\n")
		
	cat( package_version_date("sirt") , "\n" )
	cat( package_version_date("lavaan") , "\n" )	
	cat( package_version_date("lavaan.survey") , "\n" )			
	
	cat(paste0("\nFunction 'sirt::lsem.estimate', type='" , object$type,"'") , "\n\n")
		
	cat("Call:\n", paste(deparse(object$CALL), sep = "\n", collapse = "\n"), 
				"\n\n", sep = "")		
	
	cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
	cat("Computation Time:" , print(object$s2 - object$s1), "\n\n")
    
	
	cat( paste0( "Number of observations = " , round(object$N,digits) ) , "\n")
	if ( object$type == "LSEM"){
		cat( paste0( "Bandwidth factor = " , round(object$h,digits) ) , "\n")
		cat( paste0( "Bandwidth = " , round(object$bw,digits) ) , "\n")
		cat( paste0( "Number of focal points for moderator = " , 
							length(object$moderator.grid ) ) , "\n")
								}

	if ( object$type == "MGM"){
#		cat( paste0( "Bandwidth factor = " , round(object$h,digits) ) , "\n")
#		cat( paste0( "Bandwidth = " , round(object$bw,digits) ) , "\n")
		cat( paste0( "Number of groups for moderator = " , 
							length(object$moderator.grid ) ) , "\n")
								}								
								
	cat("\nlavaan Model\n")
	cat(object$lavmodel)						
						
    cat("\n\n")
	cat("Parameter Estimate Summary\n\n")
	obji <- object$parameters_summary
	VV <- ncol(obji)
	for (vv in 2:VV){
		obji[,vv] <- round( obji[,vv] , digits )
					}
	print(obji)
	
    cat("\n")
	cat("Distribution of Moderator: Density and Effective Sample Size\n\n")
	
	obji <- object$moderator.density
	VV <- ncol(obji)
	for (vv in 1:VV){
		obji[,vv] <- round( obji[,vv] , digits )
					}
	print(obji)		
	cat("\n")
	obji <- object$moderator.stat	
	VV <- ncol(obji)
	for (vv in 2:VV){
		obji[,vv] <- round( obji[,vv] , digits )
					}
	print(obji)	
	
	# close file
	CDM::csink(file)
	
	}
#############################################
