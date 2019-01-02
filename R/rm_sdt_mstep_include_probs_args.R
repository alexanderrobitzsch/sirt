## File Name: rm_sdt_mstep_include_probs_args.R
## File Version: 0.04

rm_sdt_mstep_include_probs_args <- function(probs_args, parm_list,
    update_probs_args )
{
    for ( uu in update_probs_args){
        probs_args[[ uu ]] <- parm_list[[ uu ]]
    }
    return(probs_args)
}
