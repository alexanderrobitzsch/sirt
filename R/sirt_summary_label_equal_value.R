## File Name: sirt_summary_label_equal_value.R
## File Version: 0.06

sirt_summary_label_equal_value <- function( label, value, label2="", digits=NULL )
{
    if ( ! is.null(digits) ){
        value <- round( value, digits)
    }
    res <- paste0( label, " ", "=", " ", value, " ", label2 )
    return(res)
}
