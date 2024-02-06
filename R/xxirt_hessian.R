## File Name: xxirt_hessian.R
## File Version: 0.579

#--- computation of hessian matrix
xxirt_hessian <- function( object, h=1e-4, use_shortcut=TRUE )
{
    item_list <- object$item_list
    items <- object$items
    Theta <- object$Theta
    ncat <- object$ncat
    partable <- object$partable
    partable_index <- object$partable_index
    dat <- as.matrix(object$dat)
    dat_resp <- object$dat_resp
    dat_resp_bool <- object$dat_resp_bool
    dat1 <- object$dat1
    dat1_resp <- object$dat1_resp
    resp_index <- object$resp_index
    G <- object$G
    group <- object$group
    maxK <- object$maxK
    group_index <- object$group_index
    weights <- object$weights
    customTheta <- object$customTheta
    par_items <- object$par_items
    par_Theta <- object$par_Theta
    NPI <- length(par_items)
    NPT <- length(par_Theta)

    #** preliminary computations
    probs_items <- xxirt_compute_itemprobs( item_list=item_list,
                                items=items, Theta=Theta, ncat=ncat,
                                partable=partable, partable_index=partable_index )

    par1 <- xxirt_partable_extract_freeParameters( partable=partable )

    p.xi.aj <- xxirt_compute_likelihood( probs_items=probs_items, dat=dat,
                             resp_index=resp_index, dat_resp_bool=dat_resp_bool )
    prior_Theta <- xxirt_compute_priorDistribution( Theta=Theta,
                                  customTheta=customTheta, G=G )
    prior1 <- t( prior_Theta[, group ] )
    probs_items0 <- probs_items
    p_xi_aj <- p.xi.aj

    #--- compute Hessian matrix
    par1 <- xxirt_partable_extract_freeParameters( partable=partable )
    par2 <- xxirt_parTheta_extract_freeParameters( customTheta=customTheta )
    par0 <- par <- c(par1, par2)

    #** detect whether there are item-wise parameters
    a1 <- aggregate( partable$itemnr, list(partable$parindex), min )
    a2 <- aggregate( partable$itemnr, list(partable$parindex), max )
    item_wise <- sum(a2[,2]> a1[,2])==0

    #***********************
    fct_irt <- function(x, prior1, par0, item_wise=FALSE, probs_items0,
                    p_xi_aj)
    {
        eps <- 1e-10
        I <- dim(probs_items0)[1]
        #*** include free paramaters in partable
        partable <- xxirt_partable_include_freeParameters( partable, x=x[ 1:NPI ] )
        #**** include parameter in customTheta
        customTheta$par[ customTheta$est ] <- x[ (NPI+1):(NPI+NPT) ]
        #*** different parameters
        ind <- which(abs(x-par0)>eps)
        max_ind <- 0
        min_ind <- 0
        IA <- 100
        if (length(ind)>0){
            max_ind <- max(ind)
            min_ind <- min(ind)
            items_active <- unique(partable[ partable$parindex %in% ind, 'itemnr' ])
            if (length(items_active)>0){
                IA <- length(items_active)
            }
        }

        #*** compute individual likelihood
        # if ( (! item_wise) | (min_ind==0) | (max_ind>=NPI+1 ) ){

        use_shortcut <- use_shortcut & (IA <=2)

        if (! use_shortcut){
            if (min_ind <=NPI){
                probs_items <- xxirt_compute_itemprobs( item_list=item_list,
                                    items=items, Theta=Theta, ncat=ncat,
                                    partable=partable, partable_index=partable_index )
            }
            p.xi.aj <- xxirt_compute_likelihood( probs_items=probs_items, dat=dat,
                             resp_index=NULL, dat_resp_bool=dat_resp_bool )
        } else {
            itemnr <- partable[ partable$parindex==min_ind, 'itemnr' ]
            itemnr2 <- partable[ partable$parindex==max_ind, 'itemnr' ]
            maxK <- dim(probs_items)[2]
            TP <- dim(probs_items)[3]
            eps <- 1e-16
            use_itemnr2 <- FALSE
            if (length(itemnr2)>0){
                if (itemnr2!=itemnr){
                    use_itemnr2 <- TRUE
                }
            } else {
                if (itemnr>1){
                    itemnr2 <- 1
                } else {
                    itemnr2 <- 2
                }
            }

            item_index <- c(itemnr, itemnr2)
            if (item_index[2]==item_index[1]){
                if (item_index[1]==1){
                    item_index[2] <- 2
                } else {
                    item_index[2] <- 1
                }
            }
            probs_items_temp <- xxirt_compute_itemprobs( item_list=item_list,
                                    items=items, Theta=Theta, ncat=ncat,
                                    partable=partable, partable_index=partable_index,
                                    item_index=item_index)
            probs_ratio1 <- probs_items0

            p0 <- probs_items0[item_index,,]+eps
            probs_ratio1[item_index,,] <- probs_items_temp / p0
            probs_ratio1 <- matrix( probs_ratio1, nrow=I, ncol=maxK*TP )

            p.xi.aj <- sirt_rcpp_xxirt_hessian_reduced_probs(dat=dat,
                            dat_resp_bool=dat_resp_bool, probs_ratio=probs_ratio1,
                            TP=TP, maxK=maxK, itemnr=itemnr-1, itemnr2=itemnr2-1,
                            use_itemnr2=use_itemnr2, p_xi_aj=p_xi_aj)

        }

        #*** compute prior distribution
        if (max_ind >=NPI+1){
            prior_Theta <- xxirt_compute_priorDistribution( Theta=Theta,
                                  customTheta=customTheta, G=G )
            prior1 <- t( prior_Theta[, group ] )
        }
        dev <- xxirt_hessian_compute_loglike(p.xi.aj=p.xi.aj, prior1=prior1,
                    weights=weights)
        return(dev)
    }
    #*******************************

    hess <- CDM::numerical_Hessian( par=par, FUN=fct_irt, h=h, prior1=prior1,
                        par0=par, item_wise=item_wise, probs_items0=probs_items0,
                        p_xi_aj=p_xi_aj)
    rownames(hess) <- colnames(hess) <- names(par)
    return(hess)
}
