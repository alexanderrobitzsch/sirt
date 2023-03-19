## File Name: linking_haebara_summary_optimization.R
## File Version: 0.091

linking_haebara_summary_optimization <- function(object, digits)
{
    cat('Distance function type', '=', object$dist, '\n')
    if (object$dist=='L1'){
        cat('Epsilon Value', '=', object$eps, '\n')
    }
    cat('Optimization Function Value', '=', round(object$res_optim$value, digits), '\n')
    cat('Optimizer', '=', object$res_optim$optimizer, '\n')
    cat('use_rcpp', '=', object$use_rcpp, '\n')
    cat('Number of iterations', '=', object$res_optim$iter, '\n')
    cat('Converged', '=', object$res_optim$converged, '\n')
}
