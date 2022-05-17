## File Name: mgsem_moments_derivative_parameter.R
## File Version: 0.348

mgsem_moments_derivative_parameter <- function(est, est_add=NULL, type,
        i1, i2, h, is_B, eps=1e-12, num_approx=FALSE)
{
    calc_Sigma <- TRUE
    calc_Mu <- TRUE
    symm <- FALSE
    I <- nrow(est$LAM)
    log0 <- Sigma_der_logical <- matrix(FALSE, nrow=I, ncol=I)

    if (type %in% c("PHI","PSI")){
        symm <- TRUE
        calc_Mu <- FALSE
    }
    if (type %in% c("NU","ALPHA")){
        calc_Sigma <- FALSE
    }


    if (is_B){
        num_approx <- TRUE
    }
    if (type %in% c("LAM","PHI")){
        num_approx <- FALSE
    }

    # leave always numerical approximation because est is a function of
    # different components

    #*** numerical approximation
    #if (TRUE){
    if (num_approx){
        est1 <- est
        est1[[type]] <- mgsem_add_increment(x=est[[type]],h=h,i1=i1,i2=i2, symm=symm)
        est2 <- est
        est2[[type]] <- mgsem_add_increment(x=est[[type]],h=-h,i1=i1,i2=i2, symm=symm)

        ## here: est_group=est_between + est_within_group
        # add to est_add!

        # est3 <- mgsem_add_list_entries(list1=est1, add_list=est2, elements=type)

        #* compute implied moments
        # implied <- mgsem_compute_model_implied_moments(est=est)
        implied1 <- mgsem_compute_model_implied_moments(est=est1, is_B=is_B,
                        calc_Mu=calc_Mu, calc_Sigma=calc_Sigma)
        implied2 <- mgsem_compute_model_implied_moments(est=est2, is_B=is_B,
                        calc_Mu=calc_Mu, calc_Sigma=calc_Sigma)

        # compute numerical derivatives with respect to Mu and Sigma
        if (calc_Mu){
            Mu_der <- (implied1$Mu-implied2$Mu)/(2*h)
        } else {
            Mu_der <- matrix(0, nrow=I, ncol=1)
        }
        if (calc_Sigma){
            Sigma_der <- (implied1$Sigma-implied2$Sigma)/(2*h)
            if (type %in% c("B","PHI")){
                Sigma_der_logical <- mgsem_differ_from_zero(x=Sigma_der, eps=eps)
            }
            if (type %in% c("PSI")){
                Sigma_der_logical[i1,i2] <- TRUE
                Sigma_der_logical[i2,i1] <- TRUE
            }
            if (type %in% c("LAM")){
                Sigma_der_logical[i1,] <- TRUE
                Sigma_der_logical[,i1] <- TRUE
            }
        } else {
            Sigma_der <- matrix(0, nrow=I, ncol=I)
        }


    } else {
    ## this part must be cleaned!!
        I <- nrow(est$LAM)
        mat0 <- matrix(0, nrow=I, ncol=I)
        Mu_der <- matrix(0, nrow=I, ncol=1)
        Sigma_der <- mat0
        Sigma_der_logical <- log0

        if (type=="ALPHA"){
            Mu_der <- est$LAM[,i1,drop=FALSE]
        }
        if (type=="NU"){
            Mu_der[i1,1] <- 1
        }
        if (type=="PSI"){
            Sigma_der <- mgsem_add_increment(x=mat0,h=1,i1=i1,i2=i2, symm=symm)
            Sigma_der_logical[i1,i2] <- TRUE
            Sigma_der_logical[i2,i1] <- TRUE
        }
        if (type=="LAM"){
            Mu_der[i1,i2] <- est$ALPHA[i2,1]
            D1 <- matrix(0, nrow=I, ncol=ncol(est$LAM))
            D1[i1,i2] <- 1
            LP <- est$LAM %*% est$PHI
            # Sigma_der <- D1 %*% est$PHI %*% t(est$LAM) + est$LAM %*% est$PHI %*% t(D1)
            # Sigma_der <- D1 %*% t(LP) + LP %*% t(D1)
            # G1 <- D1 %*% t(LP)
            G1 <- tcrossprod(D1, LP)
            Sigma_der <- G1 + t(G1)
            Sigma_der_logical[i1,] <- TRUE
            Sigma_der_logical[,i1] <- TRUE
        }
        if (type=="PHI"){
            D <- nrow(est$PHI)
            D1 <- matrix(0, nrow=D, ncol=D)
            D1[i1,i2] <- 1
            D1[i2,i1] <- 1
            # Sigma_der <- est$LAM %*% D1 %*% t(est$LAM)
            Sigma_der <- est$LAM %*% tcrossprod(D1, est$LAM)
            Sigma_der_logical <- mgsem_differ_from_zero(x=Sigma_der, eps=eps)
        }
    }

    #--- output
    res <- list(Mu_der=Mu_der, Sigma_der=Sigma_der, Sigma_der_logical=Sigma_der_logical,
                    calc_Sigma=calc_Sigma, calc_Mu=calc_Mu)
    return(res)
}
