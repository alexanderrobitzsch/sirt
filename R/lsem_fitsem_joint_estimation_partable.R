## File Name: lsem_fitsem_joint_estimation_partable.R
## File Version: 0.05
## File Last Change: 2022-07-27


lsem_fitsem_joint_estimation_partable <- function(lavfit, G, par_invariant=NULL,
    par_linear=NULL, par_quadratic=NULL, pw_linear=1, pw_quadratic=1)
{
    partable <- sirt_import_lavaan_parameterTable(lavfit)
    partable$start <- partable$est
    partable_joint <- lsem_fitsem_joint_estimation_prepare_partable(partable=partable,
                                G=G, par_invariant=par_invariant, par_linear=par_linear,
                                par_quadratic=par_quadratic, pw_linear=pw_linear,
                                pw_quadratic=pw_quadratic)
    return(partable_joint)
}
