## File Name: invariance_alignment_choose_fixed.R
## File Version: 0.03
## File Last Change: 2019-08-05

invariance_alignment_choose_fixed <- function(fixed, G, Gmax=6)
{
    if (is.null(fixed)){
        if (G>Gmax){
            fixed <- FALSE
        } else {
            fixed <- TRUE
        }
    }
    return(fixed)
}
