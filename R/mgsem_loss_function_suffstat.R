## File Name: mgsem_loss_function_suffstat.R
## File Version: 0.184

mgsem_loss_function_suffstat <- function(suffstat, Mu, Sigma, p=2, eps=1e-3,
        deriv=FALSE, approx_method="lp", only_deriv=FALSE, output_all=FALSE )
{
    N <- suffstat$N
    M <- suffstat$M
    S <- suffstat$S
    if (missing(Sigma)){
        res <- Mu
        Mu <- res$Mu
        Sigma <- res$Sigma
    }
    I <- length(Mu)

    #** mean structure
    t1 <- 0
    m1 <- M-Mu
    w <- suffstat$weights_M
    if(!only_deriv){
        y <- mgsem_power_fun_differentiable_approx(x=m1, p=p, eps=eps, deriv=FALSE,
                approx_method=approx_method)
        t1 <- sum(y*w)
    }
    if (deriv){
        y <- mgsem_power_fun_differentiable_approx(x=m1, p=p, eps=eps,
                deriv=TRUE, approx_method=approx_method)
        t1_der <- -as.vector(y)*w
    }

    #** covariance structure
    t2 <- NULL
    v1 <- S-Sigma
    v1 <- mgsem_vech(x=v1)
    w <- suffstat$vech_weights_S
    if (! only_deriv ){
        y <- mgsem_power_fun_differentiable_approx(x=v1, p=p, eps=eps,
                deriv=FALSE, approx_method=approx_method)
        t2 <- sum(y*w)
    }
    if (deriv){
        y <- mgsem_power_fun_differentiable_approx(x=v1, p=p, eps=eps,
                deriv=TRUE, approx_method=approx_method)
        t2_der <- -y*w
    }

    # total loss function
    res <- -N/2*(t1+t2)

    if (deriv){
        res <- list(loss_fun=res, dermean=-N/2*t1_der, dercov=-N/2*t2_der)
    }
    #-- arrange output
    return(res)
}
