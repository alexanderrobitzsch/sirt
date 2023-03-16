## File Name: regpolca_postproc_ic.R
## File Version: 0.05


regpolca_postproc_ic <- function(ic, n_reg)
{
    ic$n_reg <- n_reg
    ic$np.items <- ic$np.items - ic$n_reg
    ic <- xxirt_ic_compute_criteria(ic=ic)
    return(ic)
}
