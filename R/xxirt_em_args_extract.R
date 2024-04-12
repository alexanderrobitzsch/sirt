## File Name: xxirt_em_args_extract.R
## File Version: 0.02

xxirt_em_args_extract <- function(em_args, em_out, object)
{
    v1 <- em_out[[object]]
    if (is.null(v1)){
        v1 <- em_args[[object]]
    }
    return(v1)
}
