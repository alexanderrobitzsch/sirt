## File Name: xxirt_simulate.R
## File Version: 0.142

xxirt_simulate <- function(partable, customItems, Theta, customTheta, N=1e4,
                    method="random")
{

    #---- simulate theta distribution
    P_Theta <- customTheta$P
    args_Theta <- list(par=customTheta$par, Theta=Theta, G=1)
    probs_Theta <- do.call(what=customTheta$P, args=args_Theta)[,1]
    TP <- nrow(Theta)
    indices <- as.data.frame(matrix(NA, nrow=TP, ncol=3))
    colnames(indices) <- c('start', 'end', 'N')
    if (method=='random'){
        Theta_sim <- stats::rmultinom(n=N, size=1, prob=probs_Theta)
        N_Theta <- rowSums(Theta_sim)
        indices$N <- N_Theta
        indices$start <- c(1,cumsum(N_Theta)[-TP]+1)
        indices$end <- cumsum(N_Theta)
    }
    if (method=='quasiexact'){
        N_Theta <- N*probs_Theta
        N_Theta1 <- cumsum(N_Theta)
        N_Theta2 <- floor(N_Theta1)
        dec <- N_Theta1 - N_Theta2
        indices$end <- round(N_Theta1)
        indices$start <- c(1, indices$end[-TP]+1 )
        indices$N <- indices$end - indices$start + 1
    }
    #-- define index table for theta
    if( indices$start[2]==0){
        indices$start[1] <- 0
    }
    indices$start <- ifelse( indices$start==0 & indices$N > 0, 1, indices$start )
    indices$N <- ifelse( indices$N<1, 0, indices$N )

    #---- simulate item responses
    I <- max(partable$itemnr)
    items <- unique(paste(partable$item))
    ncat <- aggregate(partable$ncat, list(partable$itemnr), max)[,2]

    dat <- as.data.frame( matrix(NA, nrow=N, ncol=I) )
    colnames(dat) <- items

    #- simulate random numbers
    mu <- rep(0,I)
    Sigma <- diag(I)
    exact <- TRUE
    if (method=='random'){ exact <- FALSE }
    rn <- rmvn(N=N, mu=mu, Sigma=Sigma, exact=exact)
    colnames(rn) <- items
    rn <- stats::pnorm(q=rn)
    # rn <- matrix( stats::runif(N*I), nrow=N, ncol=I)

    #--- simulate items
    TP_list <- which( indices$N > 0 )

    for (ii in 1:I){

        item_ii <- items[ii]
        ncat_ii <- ncat[ii]
        partable_ii <- partable[ partable$itemnr==ii, ]
        type_ii <- paste(partable_ii$type[1])
        LC <- length(customItems)
        ll0 <- NA
        for (ll in 1:LC){
            if (customItems[[ll]]$name==type_ii){
                ll0 <- ll
                customItems_ll <- customItems[[ll]]
            }
        }
        par <- partable_ii$value
        names(par) <- paste(partable_ii$parname)
        args_fun <- list( par=par, Theta=Theta, ncat=ncat_ii)
        probs <- do.call(what=customItems_ll$P, args=args_fun)

        for (tt in TP_list){
            indices_tt <- indices[tt,]
            N_tt <- indices_tt$N
            probs_tt <- probs[tt,]
            index_tt <- seq( indices_tt$start, indices_tt$end )
            sizes_tt <- probs_tt * N_tt
            # csizes_tt <- round(cumsum(sizes_tt))
            csizes_tt <- cumsum(sizes_tt)
            probs_tt <- csizes_tt/N_tt
            rn_tt <- as.vector(rn[ index_tt, ii ])
            eps <- 1e-7
            V <- 1-eps
            probs_tt[ probs_tt > V ] <- V
            quants <- stats::quantile( rn_tt, probs=probs_tt)
            if (method=='random'){
                quants <- cumsum(probs[tt,])
            }
            vals <- rep(0,N_tt)
            rn_ii <- rn[index_tt, ii ]
            hh <- 1
            for (hh in 1:(ncat_ii-1)){
                vals <- ifelse( rn_ii > quants[hh], vals+1, vals )
            }
            dat[index_tt,ii] <- vals
        }  # end tt
    } # end ii

    #--- output
    return(dat)
}
