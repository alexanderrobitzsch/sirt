## File Name: sirt_attach_list_elements.R
## File Version: 0.06
## File Last Change: 2018-12-30

sirt_attach_list_elements <- function(x, envir)
{
    vars <- names(x)
    for (vv in vars){
        assign( vv, x[[ vv ]], envir=envir )
    }
}
