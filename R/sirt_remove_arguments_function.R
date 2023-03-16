## File Name: sirt_remove_arguments_function.R
## File Version: 0.02
## File Last Change: 2019-05-17

sirt_remove_arguments_function <- function(fun, args)
{
    fun_formals <- formals(fun=fun)
    rem <- setdiff( names(args), names(fun_formals))
    args <- sirt_remove_list_entries(list=args, rem=rem)
    return(args)
}
