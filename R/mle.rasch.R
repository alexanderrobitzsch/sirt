## File Name: mle.rasch.R
## File Version: 1.03



#---------------------------------------------------------
# Maximum Likelihood Estimation (Rasch model)

mle.rasch <- function( dat, dat.resp=1-is.na(dat), b, theta, conv=.001,
            progress=FALSE, prior_sd=NULL)
{
    theta.change <- 1
    if ( progress){ cat("\n  MLE estimation  |" ) }
    I <- length(b)
    bM <- matrix( b, nrow=length(theta), length(b), byrow=TRUE )

    # prior for ability
    if (! is.null(prior_sd) ){
        is_prior <- TRUE
    } else {
        is_prior <- FALSE
    }

    while( max( abs( theta.change) > conv )){
        # calculate P and Q
        p.ia <- stats::plogis( theta - bM )
        q.ia <- 1 - p.ia
        # Likelihood
        l1 <- rowSums( dat.resp* ( dat - p.ia ) )
        # derivative of the objective function
        f1.obj <- rowSums( - dat.resp * p.ia * q.ia  )
        # add prior
        if (is_prior){
            l1 <- l1 - theta / prior_sd^2
            f1.obj <- f1.obj - 1 / prior_sd^2
        }

        # theta change
        theta.change <- - l1 / f1.obj
        theta <- theta + theta.change
        if ( progress){  cat("-") }
    }
    res <- list( "theta"=theta, "p.ia"=p.ia )
    return(res)
}

