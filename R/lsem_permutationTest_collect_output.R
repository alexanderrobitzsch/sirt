## File Name: lsem_permutationTest_collect_output.R
## File Version: 0.02

lsem_permutationTest_collect_output <- function(res0)
{
    res0_out <- list(est=res0$parameters$est, M=res0$parameters_summary$M,
                            SD=res0$parameters_summary$SD,
                            MAD=res0$parameters_summary$MAD,
                            lin_slo=res0$parameters_summary$lin_slo)
    return(res0_out)
}
