## File Name: mcmc_plot.R
## File Version: 0.19


######################################################
# mcmc plot
mcmc_plot <- function(mcmcobj, ...)
{
    mcmcobj <- mcmc_extract_samples_first_chain(mcmcobj=mcmcobj)
    x <- list( "mcmcobj"=mcmcobj )
    x$amh_summary <- mcmc_summary(mcmcobj)
    class(x) <- "amh"
    amh_plot(x, ... )
}
