## File Name: linking_2groups_stocking_lord_grad.R
## File Version: 0.03


linking_2groups_stocking_lord_grad <- function(x, pars, Theta, wgt, type="asymm",
        pow=2, eps=1e-3, simultan=FALSE)
{
    res <- linking_2groups_stocking_lord_fun(x=x, pars=pars, Theta=Theta, wgt=wgt,
                type=type, pow=pow, eps=eps, simultan=simultan, deriv=TRUE )
    return(res)
}
