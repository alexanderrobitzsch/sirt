## File Name: lsdm_est_logist_quant.R
## File Version: 0.193


#--- Function for calculating logistic functions and probability quantiles
lsdm_est_logist_quant <- function( probcurves, theta, quantiles, wgt_theta,
        est.icc=TRUE, b=NULL, a=NULL)
{
    # estimate parameters of attribute response probcurves
    I <- nrow(probcurves)
    b0 <- NULL
    if (est.icc){
        pars.probcurves <- matrix( 0, nrow=I, ncol=5 )
        colnames(pars.probcurves) <- c('b.2PL', 'a.2PL', 'sigma.2PL', 'b.1PL',
                                            'sigma.1PL')
        rownames(pars.probcurves) <- rownames(probcurves)
        for (kk in 1:I){
            if (!is.null(b)){
                b0 <- b[kk]
                a0 <- a[kk]
            }
            pars.probcurves[kk,1:3] <- lsdm_est_logist_2pl( y=probcurves[kk,],
                                            theta=theta, wgt_theta=wgt_theta,
                                            b0=b0, a0=a0 )
            pars.probcurves[kk,4:5] <- lsdm_est_logist_rasch( y=probcurves[kk,],
                                            theta=theta, wgt_theta=wgt_theta )
        }
    }
    # quantiles of Item Response Curves (Logistic Functions)
    probcurves.quant <- sapply( quantiles, FUN=function(ql){
            sapply( 1:I, FUN=function(kk){
                    lsdm_extract_probquantile(vec=probcurves[kk,], theta=theta, quant=ql)
                } )
            } )
    probcurves.quant <- as.data.frame(probcurves.quant)
    colnames(probcurves.quant) <- paste( 'Q', 100*quantiles, sep='')
    rownames(probcurves.quant) <- rownames(probcurves)
    if (est.icc){
        pars.probcurves <- cbind( probcurves.quant, pars.probcurves )
    } else {
        pars.probcurves <- probcurves.quant
    }
    for (vv in 1:(length(quantiles))){
        pars.probcurves[,vv] <- as.numeric( pars.probcurves[,vv] )
    }
    return(pars.probcurves)
}
