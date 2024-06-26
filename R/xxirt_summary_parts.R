## File Name: xxirt_summary_parts.R
## File Version: 0.088
#
xxirt_summary_parts <- function(object, digits, len_disp=66)
{
    sirt_display_function(length=len_disp)
    #- package and R session
    sirt_summary_print_package_rsession(pack='sirt')

    #- print call
    sirt_summary_print_call(CALL=object$CALL)

    #-- print computation time
    sirt_summary_print_computation_time_s1(object=object)

    #    modeltype <- object$irtmodel
    cat( '   ', object$ic$n, 'Cases, ', object$ic$I, 'Items, ',
                    object$G, 'Group(s)',     '\n')

    sirt_display_function(length=len_disp)

    est_ml <- object$estimator=='ML'

    cat( 'Estimator','=', object$estimator, '\n' )
    if (est_ml){
        cat( 'Number of EM iterations','=', object$iter_em, '\n' )
    }
    cat( 'Number of Newton-Raphson iterations','=', object$iter_nr, '\n' )
    if (est_ml){
        cat( 'Deviance','=', round( object$deviance, 2 ), ' | ' )
        cat( 'Log Likelihood','=', round( -object$deviance/2, 2 ), '\n' )
        cat( 'Penalty function','=', round( object$pen_val, 4 ), '\n' )
    } else {
        cat( 'Optimization function','=', round( object$res_opt_nr$objective, 2), '\n')
    }
    cat( 'Number of persons','=', object$ic$n, '\n' )

    cat( 'Number of estimated parameters','=', object$ic$np, '\n' )
    cat( '  Number of estimated item parameters','=', object$ic$np.item, '\n' )
    cat( '  Number of estimated distribution parameters','=', object$ic$np.Theta, '\n\n' )

    #--- information criteria
    if (est_ml){
        rm_summary_information_criteria(object=object)
    }

    #- trait parameters
    sirt_display_function(length=len_disp)
    cat('Trait Parameters\n')
    obji <- object$customTheta$par
    sirt_summary_print_objects(obji=obji, digits=digits, from=1)

    #- item parameters
    sirt_display_function(length=len_disp)
    cat('Item Parameters \n')
    obji <- object$par_items_summary
    sirt_summary_print_objects(obji=obji, digits=digits, from=2)
}
