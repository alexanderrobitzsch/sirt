## File Name: summary.lsdm.R
## File Version: 0.199



#*** summary for LSDM function
summary.lsdm <- function( object, file=NULL, digits=3, ... )
{
    # open sink for a file
    sirt_osink(file=file)

    lsdmobj <- object
    # generate sequence for display
    display <- paste0( paste( rep("-", each=65 ), collapse="" ), "\n")
    cat(display)

    #- package and R session
    sirt_summary_print_package_rsession(pack="sirt")

    #- print call
    sirt_summary_print_call(CALL=object$CALL)

    #-- print computation time
    sirt_summary_print_computation_time_s1(object=object$time)

    cat("\nLSDM -- Least Squares Distance Method of Cognitive Validation \n")
    cat("Reference: Dimitrov, D. (2007). Applied Psychological Measurement, 31, 367-387.\n")

    cat(display)
    cat("Model Fit\n\n")
    cat("Distance function", "=", object$distance, "\n" )
    cat(paste( "Model Fit LSDM   -  Mean MAD:",
                    formatC( lsdmobj$mean.mad.lsdm0, digits=digits, width=digits+3),
                    "    Median MAD:", formatC( median(lsdmobj$item$mad.lsdm), digits=digits, width=digits+3), "\n") )
    cat(paste( "Model Fit LLTM   -  Mean MAD:", formatC( lsdmobj$mean.mad.lltm, digits=digits, width=digits+3),
                    "    Median MAD:", formatC( median(lsdmobj$item$mad.lltm), digits=digits, width=digits+3),
                    "   R^2=", format( summary(lsdmobj$lltm)$r.squared, digits=digits),   "\n") )

    cat(display)
    cat("Attribute Parameters\n\n")
    dfr.a <- data.frame( "N.Items"=colSums(lsdmobj$Qmatrix), round( lsdmobj$attr.pars, digits) )
    elim.a <- union( grep("Q", colnames(dfr.a) ), grep( "sigma", colnames(dfr.a) ) )
    print( dfr.a[, -elim.a ] )

    cat(display)
    cat("Item Parameters\n\n")
    dfr.i <- round( lsdmobj$item, digits)
    displ.i <- sort( union( which( colnames(dfr.i) %in% c(  "a.2PL", "b.2PL", "N.Items", "b.1PL" ) ),
                    grep( "mad", colnames(dfr.i))) )
    print( round( dfr.i[, displ.i ], digits) )

    cat(display)
    cat("Discrimination Parameters\n\n")
    dfr.i <- round( lsdmobj$W, digits)
    print(dfr.i)

    sirt_csink(file=file)
}
