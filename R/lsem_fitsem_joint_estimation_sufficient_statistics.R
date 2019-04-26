## File Name: lsem_fitsem_joint_estimation_sufficient_statistics.R
## File Version: 0.03

lsem_fitsem_joint_estimation_sufficient_statistics <- function(partable_joint,
    is_meanstructure, sample_stats, lavaan_est_fun, ...)
{
    wmean <- sample_stats$wmean
    wcov <- sample_stats$wcov
    Nobs <- sample_stats$Nobs
    if (is_meanstructure){
        sample_mean <- wmean
    } else {
        sample_mean <- NULL
    }
    survey.fit <- lavaan_est_fun(partable_joint, sample.cov=wcov, sample.mean=sample_mean,
                        sample.nobs=Nobs, ... )
    return(survey.fit)
}
