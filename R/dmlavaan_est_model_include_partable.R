## File Name: dmlavaan_est_model_include_partable.R
## File Version: 0.06

dmlavaan_est_model_include_partable <- function(partable, se_sw)
{
    ind <- which(partable$free>0)
    partable1 <- partable[ind,]
    partable$se_sw <- 0
    partable1$se_sw <- se_sw[ partable1$pnid ]
    partable[ind,] <- partable1
    #-- output
    return(partable)
}
