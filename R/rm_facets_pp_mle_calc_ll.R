## File Name: rm_facets_pp_mle_calc_ll.R
## File Version: 0.14
## File Last Change: 2019-01-02



# calculate individual likelihood for item ii
rm_facets_pp_mle_calc_ll <- function( probs, data, ii, eps=1e-20 )
{
    N <- nrow(data)
    probs <- log(probs)
    m1 <- matrix(1:N, nrow=N, ncol=2)
    m1[,2] <- data[,ii] + 1
    h1 <- probs[ m1 ]
    h1[ is.na(h1) ] <- 0
    return(h1)
}
