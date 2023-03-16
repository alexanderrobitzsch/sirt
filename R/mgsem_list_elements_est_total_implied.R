## File Name: mgsem_list_elements_est_total_implied.R
## File Version: 0.02
## File Last Change: 2022-02-25

mgsem_list_elements_est_total_implied <- function(model, is_B)
{
    G <- length(model)-1
    est_total0 <- list()
    implied0 <- list()
    for (gg in 1:G){
        est0 <- model[[1]]$est
        est_gg <- model[[gg+1]]$est
        est_total0[[gg]] <- mgsem_add_list_entries(list1=est0, add_list=est_gg,
                                output_list=est0)
        implied0[[gg]] <- mgsem_compute_model_implied_moments(est=est_total0[[gg]],
                                    is_B=is_B, calc_Sigma=TRUE, calc_Mu=TRUE)
    }
    #--- output
    res <- list(est_total0=est_total0, implied0=implied0)
    return(res)
}
