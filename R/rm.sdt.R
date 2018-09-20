## File Name: rm.sdt.R
## File Version: 8.875

#################################################################
# Hierarchical rater model
# MML estimation
rm.sdt <- function( dat, pid, rater,Qmatrix=NULL, theta.k=seq(-9,9,len=30),
    est.a.item=FALSE, est.c.rater="n",
    est.d.rater="n", est.mean=FALSE, est.sigma=TRUE, skillspace="normal",
    tau.item.fixed=NULL, a.item.fixed=NULL,
    d.min=.5, d.max=100,  d.start=3,  c.start=NULL, tau.start=NULL, sd.start=1,
    d.prior=c(3,100), c.prior=c(3,100), tau.prior=c(0,1000), a.prior=c(1,100),
    link_item="GPCM", max.increment=1, numdiff.parm=.00001, maxdevchange=.10,
    globconv=.001, maxiter=1000, msteps=4, mstepconv=.001, optimizer="nlminb" )
{
    #..........................................................
    CALL <- match.call()
    s1 <- Sys.time()
    theta.k0 <- theta.k
    pi.k <- sirt_dnorm_discrete(x=theta.k, mean=0, sd=sd.start)
    max.b.increment <- max.increment
    a_center_type <- 2
    PEM <- FALSE

    #-- process data
    procdata <- res <- rm_proc_data( dat=dat, rater=rater, pid=pid )
    dat2 <- as.matrix(res$dat2)
    dat2.resp <- as.matrix(res$dat2.resp)
    rater.index1 <- res$rater.index
    dataproc.vars <- res$dataproc.vars
    VV <- res$VV
    RR <- res$RR
    item.index <- res$dataproc.vars$item.index
    rater.index <- res$dataproc.vars$rater.index
    dat2.ind.resp <- res$dat2.ind.resp
    I0 <- res$I0

    #* processing for Rcpp code
    item_index <- item.index - 1
    dat2_resp <- dat2.resp==1

    deviance.history <- rep(NA, maxiter )

    # maximum categories
    maxK <- sirt_colMaxs(dat)
    K <- max( maxK )
    K1 <- K + 1
    if ( is.null(Qmatrix) ){
        Qmatrix <- matrix( 1:K, nrow=VV, ncol=K, byrow=TRUE)
    }
    TP <- length(theta.k)
    I <- VV*RR

    # define constraints on tau.item parameters
    # if not all categories are observed
    if ( is.null( tau.item.fixed )){
        tau.item.fixed <-  rm_determine_fixed_tau_parameters( K=K, maxK=maxK, VV=VV )
    }

    # starting values for item difficulties
    tau.item <- matrix( 0, nrow=VV, ncol=K )
    rownames(tau.item) <- colnames(dat)

    tau.item <- matrix( seq( -2, 2, len=K ), nrow=VV, ncol=K, byrow=TRUE )
    if ( ! is.null(tau.start) ){
        tau.item <- tau.start
    }
    if ( ! is.null(tau.item.fixed) ){
        tau.item[ tau.item.fixed[,1:2,drop=FALSE] ] <- tau.item.fixed[,3]
    }

    a.item <- rep( 1, VV )
    if (skillspace=="discrete" ){
        est.mean <- TRUE
    }
    if ( ! is.null( a.item.fixed ) ){
        est.a.item <- TRUE
        a.item[ a.item.fixed[,1] ] <- a.item.fixed[,2]
    }

    # rater parameter
    d.rater <- matrix( d.start, nrow=I, ncol=1 )
    if (is.null(c.start)){
        c.rater <- matrix( d.start*((1:K) - .5 ), nrow=I, ncol=K, byrow=TRUE )
    }
    if ( ! is.null(c.start) ){
        c.rater <- c.start
    }

    # set c.rater for fixed items to 99
    c.rater.fixed <- NULL
    if ( ! is.null( tau.item.fixed ) ){
        tau1 <- tau.item.fixed[ tau.item.fixed[,3]==99,, drop=FALSE]
        ind <- match( item.index, tau1[,1] )
        c.rater.fixed <- tau1[ ind, ]
        c.rater.fixed[,1] <- seq( 1, nrow(c.rater.fixed) )
        c.rater.fixed[,3] <- 999
        c.rater.fixed <- c.rater.fixed[ ! is.na( c.rater.fixed[,2] ), ]
        c.rater[ c.rater.fixed[,1:2] ] <- c.rater.fixed[,3]
    }

    #--- indices for derivatives
    diffindex <- rm_sdt_prepare_diffindex( item.index=item.index,
                        rater.index=rater.index, I=I, est.c.rater=est.c.rater,
                        est.d.rater=est.d.rater )

    # init standard errors
    se.d.rater <- NA*d.rater
    se.c.rater <- NA*c.rater
    se.a.item <- NA*a.item

    d.rater.incr <- 2
    tau.item.incr  <- max.b.increment
    c.rater.incr <- max.b.increment
    a.item.incr <- max.b.increment

    #*** create parameter table
    res <- rm_sdt_create_partable( item.index=item.index,
                rater.index=rater.index, est.c.rater=est.c.rater, est.d.rater=est.d.rater,
                tau.item=tau.item, c.rater=c.rater, diffindex=diffindex, tau.prior=tau.prior,
                a.prior=a.prior, d.prior=d.prior, c.prior=c.prior, est.a.item=est.a.item,
                tau.item.fixed=tau.item.fixed, a.item.fixed=a.item.fixed,
                c.rater.fixed=c.rater.fixed)
    partable_item <- res$partable_item
    partable_rater <- res$partable_rater
    par_index <- res$par_index
    pargroup_item <- res$pargroup_item
    pargroup_rater <- res$pargroup_rater

    #*** fill parameter table with initial values
    res <- rm_sdt_fill_init_partables( partable_item=partable_item, partable_rater=partable_rater,
                par_index=par_index, tau.item=tau.item, a.item=a.item, c.rater=c.rater,
                d.rater=d.rater)
    partable_item <- res$partable_item
    partable_rater <- res$partable_rater

    #*** optimization
    if (optimizer=="optim"){
        ctl <- list(maxit=msteps)
    }
    if (optimizer=="nlminb"){
        ctl <- list(iter.max=msteps)
    }

    # inits
    iter <- 0
    dev0 <- dev <- 1E300
    parm_minimal <- list(dev=dev)
    conv <- devchange <- 1000
    mu <- 0
    sigma <- 1
    disp <- "...........................................................\n"
    prob.item <- NULL
    prob.rater <- NULL

    #****************************************************
    # start EM algorithm
    while( ( ( maxdevchange < devchange ) | (globconv < conv) ) &
            ( iter < maxiter ) ){
        cat(disp)
        cat("Iteration", iter+1, "   ", paste( Sys.time() ), "\n" )

        # previous values
        d.rater0 <- d.rater
        tau.item0 <- tau.item
        dev0 <- dev
        mu0 <- mu
        sigma0 <- sigma
        a.item0 <- a.item
        c.rater0 <- c.rater

        #** probabilities item level
        if (link_item=="GPCM" ){
            prob.item <- rm_sdt_calc_probs_gpcm_rcpp( a.item=a.item, tau.item=tau.item, Qmatrix=Qmatrix,
                            theta.k=theta.k, VV=VV, K=K, TP=TP, eps=0, use_log=FALSE, as_vector=FALSE )
        }
        if (link_item=="GRM" ){
            prob.item <- rm_sdt_calc_probs_grm_item_rcpp( tau.item=tau.item,
                        a.item=a.item, theta.k=theta.k, VV=VV, K=K, TP=TP, eps=0, use_log=FALSE )
        }
        prob_item <- as.vector(prob.item)

        #** probabilities rater level
        prob.rater <- rm_sdt_calc_probs_grm_rcpp( c.rater=c.rater, d.rater=d.rater, I=I, K=K,
                                eps=0, use_log=FALSE )
        prob_rater <- as.vector(prob.rater)

        #-- calculate posterior
        res <- sirt_rcpp_rm_sdt_posterior( prob_item=prob_item, prob_rater=prob_rater, I0=I0, K1=K1,
                    TP=TP, item_index=item_index, dat2=dat2, dat2_resp=dat2_resp, I=I, pi_k0=pi.k )
        nik_item <- res$nik_item
        nik_rater <- res$nik_rater
        f.qk.yi <- res$fqkyi
        f.yi.qk <- res$fyiqk
        pi.k <- res$pi_k
        like <- res$like
        ll <- res$ll

        #-- item level
        # (item, true_category, node)
        nik.item <- array(nik_item, dim=c(I0, K1, TP) )

        #-- rater level
        # (item_rater, obs_category, true_category)
        nik.rater <- array(nik_rater, dim=c(I, K1, K1) )

        #---- item level: optimization and gradient function
        if (link_item=="GPCM" ){
            expected_loglike_item <- function(x, ...){
                post <- rm_sdt_mstep_item_function_value( x=x, par_index=par_index,
                            partable_item=partable_item, nik.item=nik.item, Qmatrix=Qmatrix,
                            theta.k=theta.k, VV=VV, K=K, TP=TP, eps=1E-10 )
                return(post)
            }
            gradient_loglike_item <- function(x, ...){
                grad_post <- rm_sdt_mstep_item_function_gradient( x=x, par_index=par_index,
                                    partable_item=partable_item, Qmatrix=Qmatrix, theta.k=theta.k, VV=VV,
                                    K=K, TP=TP, pargroup_item=pargroup_item, nik.item=nik.item,
                                    numdiff.parm=numdiff.parm, eps=1E-10 )
                return(grad_post)
            }
        }
        if (link_item=="GRM" ){
            #- arguments
            probs_fun <- rm_sdt_calc_probs_grm_item_rcpp
            K1 <- K+1
            probs_dim <- c(VV, K1, TP)
            probs_args <- list( VV=VV, K=K, theta.k=theta.k, TP=TP, eps=1E-10, use_log=TRUE )
            update_probs_args <- c('tau.item', 'a.item')
            #- optimization function
            expected_loglike_item <- function(x, ...){
                post <- rm_sdt_mstep_type_function_value( x=x, par_index=par_index,
                            partable=partable_item, type='item', probs_args=probs_args,
                            probs_fun=probs_fun, probs_dim=probs_dim,
                            update_probs_args=update_probs_args, nik=nik_item )
                return(post)
            }
            #- gradient
            gradient_loglike_item <- function(x, ...){
                grad_post <- rm_sdt_mstep_type_function_gradient( x=x,
                        par_index=par_index, partable=partable_item, type='item',
                        pargroup_type=pargroup_item, probs_args=probs_args, probs_fun=probs_fun,
                        probs_dim=probs_dim, update_probs_args=update_probs_args, nik=nik.item,
                        numdiff.parm=numdiff.parm )
                return(grad_post)
            }
        }

        #*** optimization item level
        parm0 <- rm_sdt_extract_par_from_partable(partable=partable_item)
        res <- sirt_optimizer(optimizer=optimizer, par=parm0, fn=expected_loglike_item,
                        grad=gradient_loglike_item, control=ctl)
        parm0 <- res$par
        #-- fill parameter table
        res <- rm_sdt_fill_par_to_partable( par_index=par_index, partable=partable_item,
                        parm0=parm0, type="item" )
        partable_item <- res$partable
        tau.item <- res$parm_list$tau.item
        a.item <- res$parm_list$a.item

        #---- rater level: optimization function
        expected_loglike_rater <- function(x, ...){
            post <- rm_sdt_mstep_rater_function_value( x=x, par_index=par_index,
                        partable_rater=partable_rater, I=I, K=K, nik_rater=nik_rater, eps=1E-10 )
            return(post)
        }
        #--- rater level: gradient
        gradient_loglike_rater <- function(x, ...){
            grad_post <- rm_sdt_mstep_rater_function_gradient( x=x,
                                par_index=par_index, partable_rater=partable_rater,
                                pargroup_rater=pargroup_rater, I=I, K=K, numdiff.parm=numdiff.parm,
                                nik.rater=nik.rater, eps=1E-10 )
            return(grad_post)
        }

        #*** optimization rater level
        parm0 <- rm_sdt_extract_par_from_partable(partable=partable_rater)
        res <- sirt_optimizer(optimizer=optimizer, par=parm0, fn=expected_loglike_rater,
                        grad=gradient_loglike_rater, control=ctl)
        parm0 <- res$par

        #-- fill parameter table
        res <- rm_sdt_fill_par_to_partable( par_index=par_index, partable=partable_rater,
                        parm0=parm0, type="rater" )
        partable_rater <- res$partable
        c.rater <- res$parm_list$c.rater
        d.rater <- res$parm_list$d.rater

        #-- update distribution
        res <- rm_smooth_distribution( theta.k=theta.k, pi.k=pi.k, est.mean=est.mean,
                    skillspace=skillspace, est.sigma=est.sigma, sigma=sigma )
        pi.k <- res$pi.k
        mu <- res$mu
        sigma <- res$sigma

        #-- save deviance values
        deviance.history[iter+1] <- dev <- -2*ll

        if (dev < parm_minimal$dev ){
            parm_minimal <- list( iter_opt=iter, dev=dev, d.rater=d.rater,
                            c.rater=c.rater, tau.item=tau.item, a.item=a.item,
                            mu=mu, sigma=sigma, pi.k=pi.k)
        }

        #-- convergence criteria
        conv <- max( abs(c.rater-c.rater0), abs( c.rater-c.rater0),
                    abs( tau.item0-tau.item), abs( a.item - a.item0 ) )
        iter <- iter+1
        devchange <- abs( ( dev - dev0 ) / dev0  )

        #-- print progress
        res <- rm_sdt_print_progress( dev=dev, dev0=dev0, c.rater=c.rater, c.rater0=c.rater0,
                    d.rater=d.rater, d.rater0=d.rater0, tau.item=tau.item, tau.item0=tau.item0,
                    a.item=a.item, a.item0=a.item0, mu=mu, sigma=sigma, iter=iter )
    }
    #---------------------------- end EM algorithm

    # *********
    # arrange OUTPUT

    #-- assign elements of parm_minimal
    iter_opt <- NULL
    envir <- environment()
    res <- sirt_attach_list_elements(x=parm_minimal, envir=envir)

    # c parameters
    if ( ! is.null( c.rater.fixed ) ){
        c.rater[ c.rater.fixed[,1:2] ] <- NA
    }

    #---
    # Information criteria
    ic <- rm_sdt_postproc_ic( dev=dev, dat2=dat2, VV=VV, RR=RR, est.mean=est.mean,
                est.sigma=est.sigma, partable_item=partable_item,
                partable_rater=partable_rater, skillspace=skillspace, theta.k=theta.k )

    #---
    # person parameters
    res <- rm_facets_postproc_person( dat2=dat2, dat2.resp=dat2.resp, procdata=procdata,
                maxK=maxK, RR=RR, theta.k=theta.k, f.qk.yi=f.qk.yi )
    person <- res$person
    EAP.rel <- res$EAP.rel

    #---
    # item
    if (!is.null(tau.item.fixed)){
        K <- max(maxK)
        I <- nrow(tau.item)
        for (ii in 1:I){
            if ( maxK[ii] < K ){
                for (kk in seq(maxK[ii]+1,K) ){
                    tau.item[ii,kk] <- NA
                }
            }
        }
        # se.tau.item[ tau.item.fixed[,1:2,drop=FALSE] ] <- NA
    }

    item <- data.frame( "item"=colnames(dat),
            "N"=colSums( 1-is.na(dat)),
            "M"=colMeans( dat, na.rm=TRUE ) )
    for (kk in 1:K){ item[, paste0("tau.Cat",kk) ] <- tau.item[,kk] }
    item$a <- a.item

    # latent mean and standard deviation
    me1 <- rep(NA,VV)
    sd1 <- rep(NA,VV)
    for (ii in 1:VV){
        pii <- prob.item[ii,,]
        qii <- matrix( c(0,Qmatrix[ii,]), nrow=K+1, ncol=ncol(pii) )
        me1[ii] <- sum( colSums( qii * pii ) * pi.k )
        sd1[ii] <- sqrt( sum( colSums( qii^2 * pii ) * pi.k ) - me1[ii]^2  )
    }
    item$latM <- me1
    item$latSD <- sd1

    cat("*********************************\n")
    cat("Item Parameters\n")
    sirt_summary_print_objects(obji=item, digits=3, from=2)

    #---
    # rater
    M1 <- colSums( dat2 ) / colSums( dat2.resp )
    N <- colSums( dat2.resp )
    rater <- data.frame( "item.rater"=colnames(dat2),
            "N"=N, "M"=M1, "d"=d.rater )
    for (zz in 1:(ncol(c.rater) ) ){
        rater[, paste0("c_",zz)] <- c.rater[,zz]
    }
    # transformed c parameters
    for (zz in 1:(ncol(c.rater) ) ){
        rater[, paste0("c_",zz,".trans")] <- c.rater[,zz] / d.rater
    }

    rater <- rater[ order( paste( rater$item.rater) ), ]
    rownames(rater) <- NULL
    rownames(item) <- NULL

    rt1 <- rater
    l1 <- paste(rt1$item.rater)
    l2 <- strsplit( l1, split="-" )
    rt1$item <- unlist( lapply( l2, FUN=function(uu){ uu[[1]] } ) )
    rt1$rater <- unlist( lapply( l2, FUN=function(uu){ uu[[2]] } ) )

    #*****
    # distribution
    skill.distribution <- data.frame("theta.k"=theta.k, "pi.k"=pi.k )

    #*****
    # labels
    dimnames(prob.item) <- list( colnames(dat), paste0("Cat", 0:K), NULL )
    dimnames(prob.rater) <- list( colnames(dat2), paste0("Cat", 0:K), NULL )

    cat("*********************************\n")
    cat("Rater Parameters\n")
    sirt_summary_print_objects(obji=rater, digits=3, from=2)

    cat("*********************************\n")
    cat("EAP Reliability","=", round(EAP.rel,3), "\n")

    s2 <- Sys.time()
    res <- list( deviance=dev, ic=ic, item=item, rater=rater, person=person,
                EAP.rel=EAP.rel, mu=mu, sigma=sigma, theta.k=theta.k, pi.k=pi.k, G=1,
                tau.item=tau.item, a.item=a.item, c.rater=c.rater, d.rater=d.rater,
                f.yi.qk=f.yi.qk, f.qk.yi=f.qk.yi, prob.item=prob.item, prob.rater=prob.rater,
                nik.item=nik.item, nik.rater=nik.rater, maxK=maxK, pi.k=pi.k,
                procdata=procdata, iter=iter, theta.k=theta.k, Qmatrix=Qmatrix, s1=s1, s2=s2,
                tau.item.fixed=tau.item.fixed, rater2=rt1, maxK=maxK, link_item=link_item,
                skill.distribution=skill.distribution, skillspace=skillspace,
                iter_opt=iter_opt, CALL=CALL, deviance.history=deviance.history )
    class(res) <- "rm.sdt"
    return(res)
}

