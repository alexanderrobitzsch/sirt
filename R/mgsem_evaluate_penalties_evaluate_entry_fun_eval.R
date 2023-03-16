## File Name: mgsem_evaluate_penalties_evaluate_entry_fun_eval.R
## File Version: 0.10
## File Last Change: 2022-02-25


mgsem_evaluate_penalties_evaluate_entry_fun_eval <- function(x, fun_eval,
        args_eval, h, deriv=FALSE)
{
    if (deriv){
        args_eval$x <- x+h
        val1 <- do.call(what=fun_eval, args=args_eval )
        args_eval$x <- x-h
        val2 <- do.call(what=fun_eval, args=args_eval )
        val <- (val1-val2)/(2*h)
    } else {
        args_eval$x <- x
        val <- do.call(what=fun_eval, args=args_eval )
    }
    return(val)
}
