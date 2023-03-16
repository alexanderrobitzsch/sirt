## File Name: mgsem_numerical_gradient.R
## File Version: 0.05
## File Last Change: 2022-02-25


mgsem_numerical_gradient <- function(par, FUN, h, symmetrize=FALSE, ...)
{
    res <- CDM::numerical_gradient(par=par, FUN=FUN, h=h, ...)
    if (symmetrize){
        res <- sirt_symmetrize(x=res)
        res <- sirt_add_names(x=res, names=names(par))
    }
    return(res)
}
