## File Name: lsem_fitsem_joint_estimation.R
## File Version: 0.184

lsem_fitsem_joint_estimation <- function(partable_joint,
    is_meanstructure, sample_stats, lavaan_est_fun, se,
    sufficient_statistics=TRUE, data_joint=NULL,
    group="group__", verbose=FALSE, G=NULL, pseudo_weights=NULL, ...)
{
    wmean <- sample_stats$wmean
    wcov <- sample_stats$wcov
    Nobs <- sample_stats$Nobs

    if (is_meanstructure){
        sample_mean <- wmean
    } else {
        sample_mean <- NULL
    }
    if (sufficient_statistics){
        survey.fit <- lavaan_est_fun(partable_joint, sample.cov=wcov,
                            sample.mean=sample_mean,
                            sample.nobs=Nobs, verbose=verbose, se=se, ... )
    } else {
        # pseudo weights
        res <- lsem_fitsem_raw_data_define_pseudo_weights(dat=data_joint,
                            pseudo_weights=pseudo_weights)
        dat1 <- res$dat
        ## this function seems to be buggy, test statistics cannot be computed
        survey.fit <- lavaan_est_fun(partable_joint, data=dat1, verbose=verbose,
                            group=group, se=se, test='standard',... )
    }
    if (verbose){
        cat('\nUse lsem.bootstrap() for computing standard errors!\n\n')
    }
    return(survey.fit)
}
