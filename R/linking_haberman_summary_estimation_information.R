## File Name: linking_haberman_summary_estimation_information.R
## File Version: 0.04
## File Last Change: 2019-01-14

linking_haberman_summary_estimation_information <- function(res_opt)
{
    cat("Estimation type", "=", res_opt$estimation,"\n")
    cat("Number of iterations", "=", res_opt$iter,"\n")
    cat("Used trimming factor ('BSQ','HUB')", "=", res_opt$cutoff,"\n")
    cat("Trimming factor estimated ('BSQ','HUB')", "=", res_opt$k_estimate,"\n")
    cat("Proportion retained observation ('LTS')", "=", res_opt$lts_prop,"\n")
}
