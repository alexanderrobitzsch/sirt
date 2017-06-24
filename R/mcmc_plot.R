

######################################################
# mcmc plot
mcmc_plot <- function(mcmcobj , ...){	
	x <- list( "mcmcobj" = mcmcobj )
	x$amh_summary <- mcmc_summary(mcmcobj)
	class(x) <- "amh"
	amh_plot(x, ... )
}
