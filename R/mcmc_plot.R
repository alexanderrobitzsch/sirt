## File Name: mcmc_plot.R
## File Version: 0.14
## File Last Change: 2017-03-03 19:05:48


######################################################
# mcmc plot
mcmc_plot <- function(mcmcobj , ...){	
	x <- list( "mcmcobj" = mcmcobj )
	x$amh_summary <- mcmc_summary(mcmcobj)
	class(x) <- "amh"
	amh_plot(x, ... )
}
