## File Name: modelfit.sirt.R
## File Version: 1.138


# model fit in sirt
modelfit.sirt <- function( object )
{
    #****
    # object of class tam.mml, tam.mml.2pl or tam.fa
    # Note that only dichotomous responses are allowed
    if (inherits(object,c("tam.mml","tam.mml.2pl")) ){
        mod <- object
        posterior <- mod$hwt
        probs <- mod$rprobs
        dat <- mod$resp
        dat[ mod$resp.ind==0 ] <- NA
    }

    #--- rasch.mml
    if (inherits(object,"rasch.mml") ){
        mod <- object
        posterior <- mod$f.qk.yi
        prob1 <- mod$pjk
        probs <- array( NA, dim=c( ncol(prob1), 2, nrow(prob1)) )
        probs[, 2, ] <- t(prob1)
        probs[, 1, ] <- 1 - t(prob1)
        dat <- mod$dat
    }

    #--- rasch.mirtlc
    if (inherits(object,"rasch.mirtlc")){
        mod <- object$estep.res
        posterior <- mod$f.qk.yi
        prob1 <- mod$pjk
        probs <- array( NA, dim=c( ncol(prob1), 2, nrow(prob1)) )
        probs[, 2, ] <- t(prob1)
        probs[, 1, ] <- 1 - t(prob1)
        dat <- object$dat
    }

    #--- rasch.pml
    if ( ! inherits(object,"rasch.pml") ){
        pmlobject <- NULL
    } else {
        data <- NULL
        posterior <- NULL
        probs <- NULL
        pmlobject <- object
    }

    #--- smirt
    if (inherits(object,"smirt") ){
        # note that for polytomous response data some adaptations are
        # necessary: see modelfit in the CDM package
        mod <- object
        probs <- mod$probs
        posterior <- mod$f.qk.yi
        dat <- mod$dat
    }

    #--- smirt
    if (inherits(object,"gom") ){
        mod <- object
        probs <- mod$probs
        posterior <- mod$f.qk.yi
        dat <- mod$dat
    }

    #--- rm.facets
#    if (inherits(object,"rm.facets") ){
#        mod <- object
#        probs <- mod$probs
#        posterior <- mod$f.qk.yi
#        dat <- mod$procdata$dat2.NA
#    }

    #--- mirt
    if (inherits(object, c("ConfirmatoryClass","ExploratoryClass","SingleGroupClass"))){
        mod <- object
        mod <- mirt.wrapper.posterior(mod)
        probs <- mod$probs
        posterior <- mod$f.qk.yi
        dat <- mod$dat
    }

    #--- R2noharm, noharm.sirt
    if ( inherits(object, c("R2noharm","noharm.sirt")) ){
        # exclusion criteria for noharm.sirt
        if ( inherits(object,"noharm.sirt")) {
            if ( object$estpars$estPsi > 0 ){
                stop("Model fit cannot be calculated because of correlated residuals")
            }
            if ( ! ( object$wgtm.default ) ){
                stop( paste0( "Model fit cannot be calculated because not all",
                                    "item pairs are used for estimation") )
            }
        }
        # evaluation of posterior
        mod <- R2noharm.EAP(noharmobj=object, theta.k=seq(-6, 6, len=15 ),
                        print.output=FALSE )
        probs <- aperm( mod$probs, c(1,3,2) )
        posterior <- mod$posterior
        dat <- object$dat
    }
    # calculate modelfit.cor
    if ( inherits(object,"rasch.pml") ){
        res <- modelfit.cor.sirt.pml( data=dat, posterior=posterior, probs=probs,
                        pmlobject=pmlobject)
    } else {
        res <- CDM::modelfit.cor2( data=dat, posterior=posterior, probs=probs )
    }

    #--- model output
    class(res) <- "modelfit.sirt"
    return(res)
}
