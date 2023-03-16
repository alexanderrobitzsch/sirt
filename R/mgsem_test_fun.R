## File Name: mgsem_test_fun.R
## File Version: 0.15
## File Last Change: 2022-07-27


mgsem_test_fun <- function(test, coef, opt_fun_args)
{
    if (test){
        requireNamespace("miceadds")
        #- function evaluation
        ll <- mgsem_opt_fun(x=coef, opt_fun_args=opt_fun_args)
        #- numerical gradient
        grad1 <- mgsem_grad_fun_numeric_approx(x=coef, opt_fun_args=opt_fun_args)
        #- analytical gradient
        args <- list(x=coef, opt_fun_args=opt_fun_args)
        grad <- do.call( what=mgsem_grad_fun, args=args)

        # dfr <- cbind( grad, grad1)

        #- print
        miceadds::Revalpr("ll")
        miceadds::Revalpr_maxabs("grad","grad1")
        stop()
    }
}
