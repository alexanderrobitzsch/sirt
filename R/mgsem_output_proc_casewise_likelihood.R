## File Name: mgsem_output_proc_casewise_likelihood.R
## File Version: 0.06
## File Last Change: 2022-02-25


mgsem_output_proc_casewise_likelihood <- function(data_proc, implied, estimator="ML")
{
    case_ll <- NULL
    if (estimator=="ML" & ( ! is.null(data_proc) ) ){
        requireNamespace("mvtnorm")
        N <- data_proc$N
        G <- data_proc$G
        data <- data_proc$data
        idgroup <- data_proc$idgroup
        case_ll <- rep(NA,N)
        for (gg in 1:G){
            implied_gg <- implied[[gg]]
            ind_gg <- which(idgroup==gg)
            Mu <- as.vector(implied_gg$Mu[,1])
            Sigma <- implied_gg$Sigma
            y <- mvtnorm::dmvnorm(x=data[ind_gg,], mean=Mu, sigma=Sigma, log=TRUE)
            case_ll[ind_gg] <- y
        }
    }
    #-- output
    return(case_ll)
}
