## File Name: rm_sdt_create_partable_include_priors.R
## File Version: 0.03

rm_sdt_create_partable_include_priors <- function(partable, type, obj)
{
    vars <- c("prior_M", "prior_SD")
    for (vv in 1:2){
        partable[ ( partable$type==type) & ( partable$est ), vars[vv] ] <- obj[vv]
    }
    return(partable)
}
