## File Name: sirt_summary_print_elapsed_time.R
## File Version: 0.04

sirt_summary_print_elapsed_time <- function(object)
{
    cat( 'Elapsed time pre-processing', '=', ' ')
    print(object$time$time_pre)
    cat( 'Elapsed time optimization', '=', ' ')
    print(object$time$time_opt)
    cat( 'Elapsed time post-processing', '=', ' ')
    print(object$time$time_post)
    cat('\n')
}
