## File Name: mgsem_eval_lp_penalty_vector.R
## File Version: 0.211


mgsem_eval_lp_penalty_vector <- function(x, fac, n, p, eps_approx, deriv, h, a=3.7,
    pen_type="lasso")
{
    # smoothic penalty
    if (pen_type=='smoothic'){
        val <- mgsem_smoothic_penalty(x=x, eps=eps_approx, deriv=deriv)
        val <- fac*val
    }

    # L0 penalty
    if (pen_type=='L0'){
        h <- 1e-4
        gamma <- 50
        val <- mgsem_L0_penalty(x=x, eps=eps_approx, gamma=gamma, deriv=deriv, h=h)
        val <- n*fac*val
    }

    # Lasso penalty
    if (pen_type=='lasso'){
        val <- mgsem_power_fun_differentiable_approx(x=x, p=p,
                            eps=eps_approx, deriv=deriv, approx_method='lp')
        val <- n*fac*val
    }

    # SCAD penalty
    if (pen_type=='scad'){
        if (deriv){
            val <- mgsem_power_fun_differentiable_approx(x=x+h, p=p,
                                eps=eps_approx, deriv=FALSE, approx_method='lp')
            val1 <- mgsem_scad_penalty(x=val, lambda=fac, a=a)
            val <- mgsem_power_fun_differentiable_approx(x=x-h, p=p,
                                eps=eps_approx, deriv=FALSE, approx_method='lp')
            val2 <- mgsem_scad_penalty(x=val, lambda=fac, a=a)
            val <- n*(val1-val2)/(2*h)
        } else {
            val <- mgsem_power_fun_differentiable_approx(x=x, p=p,
                                eps=eps_approx, deriv=FALSE, approx_method='lp')
            val <- n*mgsem_scad_penalty(x=val, lambda=fac, a=a)
        }
    }

    #--- output
    return(val)
}
