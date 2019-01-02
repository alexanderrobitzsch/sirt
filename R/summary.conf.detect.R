## File Name: summary.conf.detect.R
## File Version: 0.15

#**** Summary for conf.detect object
summary.conf.detect <- function( object, digits=3, file=NULL, ...)
{
    # open sink
    sirt_osink( file=file )

    display_string <- sirt_summary_print_display(symbol="-", len=65)
    cat(display_string)

    #- package and R session
    sirt_summary_print_package_rsession(pack="sirt")

    #- print call
    sirt_summary_print_call(CALL=object$CALL)

    itemcluster <- object$itemcluster
    IC <- length( unique(itemcluster) )

    des1 <- paste0( "Confirmatory DETECT Analysis with ", IC, " Item Clusters" )
    cat(des1,"\n")
    cat(paste("Bandwidth Scale:", object$bwscale, "\n" ) )

    cat(display_string)
    cat("Dimensionality Statistics \n")
    obji <- object$detect.summary
    sirt_summary_print_objects(obji=obji, digits=digits, from=1, rownames_null=FALSE)

    # close sink
    sirt_csink( file=file )
}

