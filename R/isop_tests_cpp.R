## File Name: isop_tests_cpp.R
## File Version: 0.04


isop_tests_cpp <- function( dat, dat.resp, weights, jackunits, JJ )
{
    res <- isop_tests_C( dat=dat, dat_resp=dat.resp, weights=weights,
                jackunits=jackunits, JJ=JJ )
    return(res)
}
