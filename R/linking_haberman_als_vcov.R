## File Name: linking_haberman_als_vcov.R
## File Version: 0.165

linking_haberman_als_vcov <- function( regr_resid, regr_wgt, transf_pars,
        selitems, estimation="OLS", vcov=TRUE, NS=NULL )
{
    transf_parsM <- matrix( transf_pars, nrow=nrow(regr_resid),
            ncol=ncol(regr_resid), byrow=TRUE )
    regr_resid <- regr_resid + transf_parsM
    regr_resid <- regr_resid[ selitems, ]
    regr_wgt <- regr_wgt[selitems,]
    if ( is.null(NS) ){
        NS <- ncol(regr_resid)
    }
    N <- nrow(regr_resid)
    data <- data.frame( y=matrix( regr_resid, ncol=1 ),
                        study=rep(1:NS, each=N) )
    data$wgt <- matrix( regr_wgt, ncol=1 )
    data <- stats::na.omit(data)
    for (ss in 2:NS){
        data[, paste0('X', ss) ] <- 1*(data$study==ss)
    }
    # form <- paste0( 'y ~ 0 + ', paste0('X', 2:NS, collapse=' + ' ) )
    if (vcov){
        form <- paste0( 'y ~ ', paste0('X', 2:NS, collapse=' + ' ) )
        mod <- stats::lm( stats::as.formula(form), data=data, weights=data$wgt )
        mod_vcov <- stats::vcov(mod)
        vcov <- mod_vcov[-1,-1]
        se <- sqrt( diag( mod_vcov) )[-1]
        if (estimation!='OLS'){
            vcov <- NA*vcov
            se <- NA*se
        }
    } else {
        vcov <- NULL
        se <- NA
    }
    #-- output
    res <- list( vcov=vcov, se=se )
    return(res)
}
