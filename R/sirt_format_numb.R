## File Name: sirt_format_numb.R
## File Version: 0.13
## File Last Change: 2018-12-30


#------ format numbers
sirt_format_numb <- function( x, digits )
{
    a1 <- round( x, digits ) + 10^{-(digits +1 ) }
    a1 <- substring( a1, 1, nchar(a1) - 1 )
    return(a1)
}

format.numb <- sirt_format_numb
