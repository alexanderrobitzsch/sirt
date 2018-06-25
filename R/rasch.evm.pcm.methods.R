## File Name: rasch.evm.pcm.methods.R
## File Version: 0.05

#*****************************************************************

#***********************
# variance matrix
vcov.rasch.evm.pcm <- function( object, ... ){
     return( object$vcov)
                }

#************************
# coefficients
coef.rasch.evm.pcm <- function( object, ... ){
     return( object$coef )
                }

#*****************************************************************
