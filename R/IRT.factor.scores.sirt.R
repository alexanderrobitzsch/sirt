## File Name: IRT.factor.scores.sirt.R
## File Version: 0.173


#--- rm.facets
IRT.factor.scores.rm.facets <- function( object, type="EAP", ... )
{
    # admissible factor score types
    x1 <- c('EAP','MLE','WLE')
    if ( ! ( type %in% x1 ) ){
        stop('Requested type is not supported!\n')
    }
    #**** EAP
    if ( type=='EAP'){
        res <- object$person
        res <- res[, c('pid', 'EAP', 'SE.EAP') ]
        attr(res,'type') <- type
        attr(res,'reliability') <- object$EAP.rel
    }
    #**** MLE or WLE
    if ( type %in% c('MLE','WLE') ){
        data <- object$procdata$dat2.NA
        a <- object$ipars.dat2$a
        b <- object$ipars.dat2$b
        theta0 <- object$person$EAP
        WLE <- if( type=='WLE'){ TRUE } else { FALSE }
        res <- rm_facets_pp_mle( data=data, a=a, b=b, theta=theta0,  WLE=WLE,
                        maxiter=30, maxincr=3, h=1e-3, convP=1e-3, maxval=9.99,
                        progress=TRUE )
        res <- data.frame(pid=object$person$pid, res )
        attr(res,'type') <- type
        attr(res,'reliability') <- mle.reliability( meas=res$est, se.meas=res$se )
    }
    return(res)
}

#--- rm.sdt
IRT.factor.scores.rm.sdt <- IRT.factor.scores.rm.facets

