## File Name: sirt_summary_print_call.R
## File Version: 0.06

sirt_summary_print_call <- function(CALL)
{
    s3 <- paste0(CALL, collapse=" ")
    if ( nchar(s3) < 1000 ){
        cat("Call:\n", paste(deparse(CALL), sep="\n", collapse="\n"),
                "\n\n", sep="")
    }
}
