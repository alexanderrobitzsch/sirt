## File Name: linking_2groups_haebara_grad.R
## File Version: 0.02


linking_2groups_haebara_grad <- function(x, pars, Theta, wgt, type="asymm",
        pow=2, eps=0.001, simultan=FALSE)
{
    res <- linking_2groups_haebara_fun(x=x, pars=pars, Theta=Theta, wgt=wgt, type=type,
        pow=pow, eps=eps, simultan=simultan, deriv=TRUE )
    return(res)
}
