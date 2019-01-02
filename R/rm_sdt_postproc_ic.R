## File Name: rm_sdt_postproc_ic.R
## File Version: 0.12

rm_sdt_postproc_ic <- function(dev, dat2, VV, RR, est.mean, est.sigma,
    partable_item, partable_rater, skillspace, theta.k)
{
    ic <- list( deviance=dev, n=nrow(dat2) )
    ic$VV <- VV
    ic$RR <- RR
    #****
    # skill parameters
    ic$np.skill <- 0
    if ( skillspace=="normal" ){
        ic$np.skill <- est.sigma + est.mean
    }
    if ( skillspace=="discrete" ){
        ic$np.skill <- length(theta.k) - 1
    }
    #*****
    # item parameters
    ic$np.item <- 0
    ic$np.item <- ic$np.item + sum( partable_item[ partable_item$type=="tau", "est"] )
    ic$np.item <- ic$np.item + sum( partable_item[ partable_item$type=="a", "est"] )

    #*****
    # rater parameters
    ic$np.rater <- 0
    ic$np.rater <- ic$np.rater + sum( partable_rater[ partable_rater$type=="c", "est"] )
    ic$np.rater <- ic$np.rater + sum( partable_rater[ partable_rater$type=="d", "est"] )
    ic$np <- ic$np.skill + ic$np.item + ic$np.rater

    #-- compute information criteria
    ic <- rm_ic_criteria(ic=ic)
    return(ic)
}

