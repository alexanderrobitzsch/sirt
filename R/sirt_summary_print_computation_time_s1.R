## File Name: sirt_summary_print_computation_time_s1.R
## File Version: 0.07
## File Last Change: 2023-03-08

sirt_summary_print_computation_time_s1 <- function(object)
{
    cat('Date of Analysis:', paste(object$s2 ), '\n' )
    cat('Computation Time:', print(object$s2 - object$s1), '\n\n')
}
