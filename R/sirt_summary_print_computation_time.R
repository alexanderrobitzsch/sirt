## File Name: sirt_summary_print_computation_time.R
## File Version: 0.06

sirt_summary_print_computation_time <- function( object )
{
    t1 <- object$time$start
    t2 <- object$time$end
    cat( 'Date of Analysis:', '\n' )
    cat( '   Start:', paste( t1 ), '\n' )
    cat( '   End  :', paste( t2 ), '\n' )
    cat('Computation time:', print( t2 - t1 ), '\n\n')
}
