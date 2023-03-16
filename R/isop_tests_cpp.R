## File Name: isop_tests_cpp.R
## File Version: 0.06
## File Last Change: 2018-12-30


isop_tests_cpp <- function( dat, dat.resp, weights, jackunits, JJ )
{
    res <- isop_tests_C( dat=dat, dat_resp=dat.resp, weights=weights,
                jackunits=jackunits, JJ=JJ )
    return(res)
}
