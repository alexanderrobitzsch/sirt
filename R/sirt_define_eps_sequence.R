## File Name: sirt_define_eps_sequence.R
## File Version: 0.01
## File Last Change: 2020-03-29

sirt_define_eps_sequence <- function(eps, eps_vec)
{
    eps_vec <- unique(c(eps_vec, eps))
    eps_vec <- eps_vec[ order(eps_vec, decreasing=TRUE) ]
    eps_vec <- eps_vec[ eps_vec >=eps ]
    return(eps_vec)
}
