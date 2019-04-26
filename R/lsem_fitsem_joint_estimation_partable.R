## File Name: lsem_fitsem_joint_estimation_partable.R
## File Version: 0.03


lsem_fitsem_joint_estimation_partable <- function(lavfit, G, par_invariant=NULL,
    par_linear=NULL, par_quadratic=NULL)
{
    partable <- sirt_import_lavaan_parameterTable(lavfit)
    partable$start <- partable$est
    partable_joint <- lsem_fitsem_joint_estimation_prepare_partable(partable=partable,
                                G=G, par_invariant=par_invariant, par_linear=par_linear,
                                par_quadratic=par_quadratic)
    return(partable_joint)
}
