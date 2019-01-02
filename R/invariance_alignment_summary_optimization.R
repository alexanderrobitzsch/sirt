## File Name: invariance_alignment_summary_optimization.R
## File Version: 0.05

invariance_alignment_summary_optimization <- function(object, digits)
{
    align.pow <- object$align.pow
    align.scale <- object$align.scale
    cat("Optimization Function Value", "=", round(object$fopt[1],digits), "\n" )
    cat("Optimizer", "=", object$res_optim$optimizer, "\n" )
    cat("Number of iterations", "=", object$res_optim$iter, "\n" )
    cat("Converged", "=", object$res_optim$converged, "\n" )
    cat("Alignment Power Values","=", round(align.pow[1], digits),
                            round(align.pow[2], digits), "\n")
    cat("Alignment Scale Values","=", round(align.scale[1], digits),
                            round(align.scale[2], digits), "\n")
    cat("Epsilon Value", "=", object$eps, "\n" )
}
