## File Name: mirt.wrapper.posterior.R
## File Version: 0.370


# calculation of posterior
mirt.wrapper.posterior <- function( mirt.obj, weights=NULL, group=NULL )
{
    TAM::require_namespace_msg("mirt")
    mirt_class <- class(mirt.obj)
    is_mg <- mirt_class=="MultipleGroupClass"
    data_group <- mirt.obj@Data$group
    orig_data <- mirt.obj@Data$data
    item_names <- colnames(orig_data)
    N_orig <- nrow(orig_data)
    G <- length(unique(data_group))
    N <- length(data_group)
    ind_group <- 1:N
    id_group <- 1
    groupNames <- NULL
    if (!is.null(group)){
        groupNames <- mirt.obj@Data$groupNames
        ind_group <- which(paste(data_group)==group)
        id_group <- which(groupNames==group)
    }
    mobj <- mirt.obj
    # extract theta
    Theta <- mobj@Model$Theta
    # extract theta distribution
    pi.k <- mobj@Internals$Prior[[id_group]]

    # extract full data
    fulldata <- mobj@Data$fulldata[[id_group]]
    #load items
    I <- ncol( mobj@Data$data )
    items <- vector('list', I)
    for(ii in 1:I){
        items[[ii]] <- mirt::extract.item(x=mobj, item=ii, group=id_group)
    }
    # check whether prodlist exists
    prodlist <- mobj@Model$prodlist
    Theta1 <- Theta
    if ( length(prodlist) > 0 ){
        Theta1 <- mirt_prodterms(theta0=Theta, prodlist=prodlist)
    }
    # item-wise probabilities for each Theta
    traces <- Probtrace_sirt(items=items, Theta=Theta1)
    # log-Likelihood
    f.yi.qk <- exp( fulldata %*% t(log(traces)) )
    # compute individual posterior
    N <- nrow(fulldata)
    TP <- length(pi.k)
    piM <- matrix( pi.k, nrow=N, ncol=TP, byrow=TRUE )
    f.qk.yi <- f.yi.qk * piM
    f.qk.yi <- f.qk.yi / matrix( rowSums(f.qk.yi), nrow=N, ncol=TP, byrow=FALSE )
    # maximum category
    maxK <- apply( mobj@Data$data, 2, max, na.rm=TRUE)+1
    resp <- mobj@Data$data
    resp <- resp[ ind_group, ]
    resp.ind <- 1 - is.na(resp)
    resp[ resp.ind==0 ] <- 0
    # calculate counts
    if (is.null(weights) ){
        pweights <- rep(1,N)
    } else {
        pweights <- weights
    }
    # Theta is only used for calculating dimension size
    n.ik <- mirt.wrapper.calc.counts( resp=resp, theta=Theta, resp.ind=resp.ind,
                group=group, maxK=max(maxK), pweights=pweights, hwt=f.qk.yi )
    dimnames(n.ik)[[2]] <- item_names
    exp_counts <- aperm(n.ik, perm=c(2,3,1,4))
    probs <- traces
    probs <- array( probs, dim=c(TP,max(maxK),I) )
    probs <- aperm( probs, perm=c(3,2,1) )
    dimnames(probs)[[1]] <- item_names
    # result list
    res <- list( theta.k=Theta, pi.k=pi.k, f.yi.qk=f.yi.qk, f.qk.yi=f.qk.yi,
            n.ik=n.ik, exp_counts=exp_counts, probs=probs, N=N, TP=TP, I=I, data=resp,
            maxK=maxK, G=G, groupNames=groupNames, group=group,
            orig_data=orig_data, ind_group=ind_group, N_orig=N_orig,
            mirt_class=mirt_class, pweights=pweights)
    class(res) <- "mirt"
    return(res)
}

