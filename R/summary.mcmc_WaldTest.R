## File Name: summary.mcmc_WaldTest.R
## File Version: 0.01
## File Last Change: 2019-04-26

# summary of Wald Test based on MCMC output
summary.mcmc_WaldTest <- function( object, digits=3, ... )
{
    cat("Wald Test\n")
    W1 <- sprintf( paste0("%.", digits, "f" ), object$chisq_stat["chi2"] )

    v1 <- paste0("Chi^2=",  W1, ", df=", object$chisq_stat["df"])
    v1 <- paste0( v1, ", p=", sprintf( paste0("%.", digits, "f" ),
                    object$chisq_stat["p"] ) )
    cat(v1)

    cat("\n\nSummary Hypotheses\n")
    obji <- object$hypotheses_summary
    vars <- c("parameter","MAP","SD", "Q2.5", "Q97.5", "Rhat","SERatio",
                    "effSize" )
    obji <- obji[,vars]
    NO <- ncol(obji)
    obji[,NO] <- round(obji[,NO])
    sirt_summary_print_objects(obji=obji, digits=digits, from=2)
}
