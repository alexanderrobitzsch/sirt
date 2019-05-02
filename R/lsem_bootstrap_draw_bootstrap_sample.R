## File Name: lsem_bootstrap_draw_bootstrap_sample.R
## File Version: 0.02

lsem_bootstrap_draw_bootstrap_sample <- function(data, sampling_weights,
    lsem_args)
{
    lsem_args1 <- lsem_args
    N <- nrow(data)
    ind <- sample(1:N, N, replace=TRUE)
    lsem_args1$data <- data[ind,]
    lsem_args1$sampling_weights <- sampling_weights[ind]
    return(lsem_args1)
}
