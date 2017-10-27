## File Name: xxirt_createThetaDistribution.R
## File Version: 0.08

######################################################################
xxirt_createThetaDistribution <- function( par , est , P , prior=NULL,
		prior_par1 = NULL , prior_par2 = NULL	)
{
    res <- list()
    res$par <- par
    res$est <- est
    res$P <- P
	NP <- length(par)
    res$prior <- prior
	res$prior_par1 <- prior_par1
	res$prior_par2 <- prior_par2
    class(res) <- "ThetaDistribution"
    return(res)
}
######################################################################				
