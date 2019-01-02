## File Name: summary.linking.robust.R
## File Version: 0.04


summary.linking.robust <- function( object, ... )
{
    kmax <- length(object$k.robust)
    cat("Robust linking with trimmed mean\n\n")
    cat(paste0("Number of items", "=", object$I, "\n" ) )
    cat(paste0( "Optimal trimming parameter k", "=", round( object$kopt, 4 ), " | "))
    cat(paste0( " non-robust parameter k", 0, " \n"))
    cat(paste0( "Linking constant=", round( object$meanpars.kopt, 4 ), " | "))
    cat(paste0( " non-robust estimate=", round( object$meanpars[ 1 ], 4 ), " \n"))
    cat(paste0( "Standard error=", round( object$se.kopt, 4 ), " | "))
    cat(paste0( " non-robust estimate=", round( object$se[1], 4 ), " \n"))
    cat(paste0( "DIF SD: MAD=", round( object$mad, 4 ), " (robust) | "))
    cat(paste0( "SD=", round( object$sd, 4 ), " (non-robust) \n"))
}
