## File Name: xxirt_createDiscItem.R
## File Version: 0.08

xxirt_createDiscItem <- function( name , par , est , P , lower=-Inf , 
            upper=Inf , prior=NULL , prior_par1=NULL ,
            prior_par2 = NULL)
{
    res <- list()
    res$name <- name
    res$par <- par
    res$est <- est
    res$P <- P
    res$lower <- lower
    res$upper <- upper
    res$prior <- prior
    res$prior_par1 <- prior_par1
    res$prior_par2 <- prior_par2
    class(res) <- "DiscItem"
    return(res)
}
