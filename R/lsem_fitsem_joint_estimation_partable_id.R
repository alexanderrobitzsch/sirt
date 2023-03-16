## File Name: lsem_fitsem_joint_estimation_partable_id.R
## File Version: 0.01
## File Last Change: 2019-04-26

lsem_fitsem_joint_estimation_partable_id <- function(partable_gg, partable_mg, vv)
{
    partable_gg[,vv] <- partable_gg[,vv] + max(partable_mg[,vv], na.rm=TRUE)
    return(partable_gg)
}
