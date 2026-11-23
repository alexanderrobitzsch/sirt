## File Name: split_syn_string_vec.R
## File Version: 0.02

#-- split_syn_string vectorized input
split_syn_string_vec <- function( syn, vecstr )
{
    for (vv in vecstr){
        syn <- split_syn_string( syn, vv )
    }
    return(syn)
}
