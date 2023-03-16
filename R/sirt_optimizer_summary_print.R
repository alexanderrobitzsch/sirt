## File Name: sirt_optimizer_summary_print.R
## File Version: 0.141

sirt_optimizer_summary_print <- function(res, msg="Information about optimization")
{
    digits <- 6
    cat('\n---', msg, '---\n\n')
    cat( 'Optimizer', '=', res$optimizer, '\n')
    cat( 'Converged', '=', res$converged, '\n')
    cat('Optimization Function Value', '=', round(res$value,digits), '\n' )
    cat( 'Number of iterations', '=', res$iter, '\n')
    cat( 'Elapsed time', '=', ' ')
    print(res$time)
    cat('\n')
}
