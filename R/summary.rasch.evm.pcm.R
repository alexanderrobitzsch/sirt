## File Name: summary.rasch.evm.pcm.R
## File Version: 0.18
## File Last Change: 2018-12-30

summary.rasch.evm.pcm <- function( object, digits=3, file=NULL, ... )
{
    # open sink for a file
    sirt_osink(file=file)

    display_string <- sirt_summary_print_display(symbol="-", len=65)
    cat(display_string)

    #- package and R session
    sirt_summary_print_package_rsession(pack="sirt")

    #- print call
    sirt_summary_print_call(CALL=object$CALL)

    #-- print computation time
    sirt_summary_print_computation_time_s1(object=object)

    cat("Partial Credit Model (estimated with eigenvector method) \n\n")
    options(scipen=999)
    cat( "Number of groups:", object$desc$G[1], "\n" )
    cat( "Number of persons:", object$desc$Nobs, "\n" )
    cat( "Weighted number of persons:", round( (object$desc$sum.weights), 2 ), "\n" )
    cat( "Number of items:", object$desc$N.items[1], "\n" )
    cat( paste0("Number of items per parameter: ",
            "M=", round( object$desc$M.Nitems[1],1),
            " | SD=", round( object$desc$SD.Nitems[1],1),
            " | Min=", object$desc$min.Nitems[1],
            " | Max=", object$desc$max.Nitems[1],
            "\n" )
                )
    cat( "Number of parameters:", sum(object$desc$N.pars), "\n\n" )

    cat(display_string)
    cat("Item Parameters \n\n")
    obji <- object$item
    sirt_summary_print_objects(obji=obji, digits=digits, from=4, rownames_null=FALSE)

    if (! is.null( object$difstats) ){
        cat(display_string)
        cat("DIF Tests \n\n")
        obji <- object$difstats
        sirt_summary_print_objects(obji=obji, digits=digits, from=2, rownames_null=FALSE)
    }
    options(scipen=0)

    sirt_csink(file=file)
}
