## File Name: invariance_alignment_choose_fixed.R
## File Version: 0.03

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
