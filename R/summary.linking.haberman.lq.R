## File Name: summary.linking.haberman.lq.R
## File Version: 0.123



summary.linking.haberman.lq <- function( object, digits=3, file=NULL, ... )
{
    # open sink
    sirt_osink(file=file)

    object$linking_slopes <- TRUE

    display_string <- sirt_summary_print_display(symbol='-', len=65)
    cat(display_string)

    #- package and R session
    sirt_summary_print_package_rsession(pack='sirt')

    cat('\n')
    cat(object$description, '\n')

    #- print call
    sirt_summary_print_call(CALL=object$CALL)

    #-- print computation time
    sirt_summary_print_computation_time_s1(object=object$time)

    cat('Converged', '=', object$converged, '\n' )
    cat('Estimated power values','=', object$est_pow, '\n')
    cat('Loss function power slopes','=', object$pow_slopes, '\n')
    cat('Loss function power intercepts','=', object$pow_intercepts, '\n')
    cat('Epsilon Value', '=', object$eps, '\n' )
    cat('a_log', '=', object$a_log, '\n' )
    cat('use_nu', '=', object$use_nu, '\n' )

    # cat(display_string)
    # cat('Transformation parameters (Haberman linking)\n')
    # obji <- object$transf.pars
    # sirt_summary_print_objects(obji=obji, digits=digits, from=2, rownames_null=FALSE)

    cat('\nLinear transformation for item parameters a and b\n')
    obji <- object$transf.itempars
    sirt_summary_print_objects(obji=obji, digits=digits, from=2, rownames_null=FALSE)

    cat('\nLinear transformation for person parameters theta\n')
    obji <- object$transf.personpars
    sirt_summary_print_objects(obji=obji, digits=digits, from=2, rownames_null=FALSE)

    if (FALSE){
        cat('\n-----------------------------------------------------------------\n')
        if ( ! object$linking_slopes ){
            cat('Estimated DIF effects in logarithms of item slopes \n')
            obji <- object$a.resid[ object$selitems, ]
            sirt_summary_print_objects(obji=obji, digits=digits, from=1,
                            rownames_null=FALSE)
            cat('\n')
        }
        cat('Estimated DIF effects of item intercepts \n')
        obji <- object$b.resid[ object$selitems, ]
        sirt_summary_print_objects(obji=obji, digits=digits, from=1,
                            rownames_null=FALSE)

        cat('\n-----------------------------------------------------------------\n')
        digits2 <- 1
        if ( ! object$linking_slopes ){
            cat('Used items in linking item slopes \n')
            obji <- object$a.item_stat
            sirt_summary_print_objects(obji=obji, digits=digits, from=1,
                            rownames_null=FALSE)
            cat('\n')
        }
        cat('Used items in linking item intercepts \n')
        obji <- object$b.item_stat
        sirt_summary_print_objects(obji=obji, digits=digits, from=1,
                            rownames_null=FALSE)
    }

    # close sink
    sirt_csink( file=file )
}
