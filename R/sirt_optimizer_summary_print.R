## File Name: sirt_optimizer_summary_print.R
## File Version: 0.09

sirt_optimizer_summary_print <- function(res)
{
    digits <- 6
    cat("\n--- Information about optimization ---\n\n")
    cat( "Optimizer", "=", res$optimizer, "\n")
    cat( "Converged", "=", res$converged, "\n")
    cat("Optimization Function Value", "=", round(res$value,digits), "\n" )
    cat( "Number of iterations", "=", res$iter, "\n")
    cat( "Elapsed time", "=", " ")
    print(res$time)
    cat("\n")
}
