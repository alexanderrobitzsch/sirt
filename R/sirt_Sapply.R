## File Name: sirt_Sapply.R
## File Version: 0.075


sirt_Sapply <- function(...)
{
    args <- list(...)
    fun <- utils::getFromNamespace(x='mySapply', ns='mirt')
    return( do.call(what=fun, args=args) )
}
