## File Name: mgsem_eval_lp_penalty_matrix.R
## File Version: 0.083

mgsem_eval_lp_penalty_matrix <- function(x, fac, p, n, h, eps_approx,
    pen_type="lasso", a_scad=3.7)
{
    x1 <- x
    I1 <- length(x1)
    y <- matrix(x1, nrow=I1, ncol=I1)-sirt_matrix2(x=x1, nrow=I1)
    y <- mgsem_power_fun_differentiable_approx(x=y, p=p,
                            eps=eps_approx, deriv=FALSE, approx_method="lp")
    if (pen_type=="lasso"){
        val <- fac*y
    }
    if (pen_type=="scad"){
        val <- mgsem_scad_penalty(x=y, lambda=fac, a=a_scad)
    }
    res <- sum(n*val)
    return(res)
}
