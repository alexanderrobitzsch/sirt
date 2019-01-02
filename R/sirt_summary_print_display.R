## File Name: sirt_summary_print_display.R
## File Version: 0.03

sirt_summary_print_display <- function(symbol="-", len=65)
{
    res <- paste0( paste0(rep(symbol, len), collapse=""), "\n" )
    return(res)
}
