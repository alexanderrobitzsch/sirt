## File Name: xxirt_prepare_response_data.R
## File Version: 0.095

xxirt_prepare_response_data <- function(G, group_index, weights, dat1,
            dat_resp, maxK )
{
    N <- nrow(dat_resp)
    I <- ncol(dat_resp)
    dat1_resp <- array(0, dim=c(N,I,maxK) )
    for (gg in 1L:G){
        ind_gg <- group_index[[gg]]
        for (kk in 1L:maxK){
            dat1_resp[ind_gg,,kk] <-  weights[ind_gg] * ( dat1[ind_gg, ]==(kk-1) ) *
                    ( dat_resp[ind_gg, ] )
        }  # end kk
    }  # end gg
    return(dat1_resp)
}
