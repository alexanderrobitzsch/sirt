## File Name: decategorize.R
## File Version: 0.04
## File Last Change: 2017-01-18 11:02:46

###########################################################
decategorize <- function( dat , categ_design=NULL ){

	# preliminarites
	dat4 <- dat3 <- dat
	dfr <- categ_design 
	
	#****************************
	# handle categories
	if ( ! is.null( dfr ) ){	
		vars <- sort( unique( paste( dfr$variable )))
		VV <- length(vars)
		for (vv in 1:VV){
			# vv <- 3
			dfr.vv <- dfr[ paste(dfr$variable) == vars[vv] , ]
			dat4[, vars[vv] ] <- dfr.vv[ match( dat3[,vars[vv]] , dfr.vv$recode ) , "orig"] 
					}				
			}	
	#***************************				
	return(dat4)
		}	
#################################################################
