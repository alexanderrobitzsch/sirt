## File Name: hh_linking_apply_precision_weights.R
## File Version: 0.02


hh_linking_apply_precision_weights <- function(V, prec_weights="none")
{
    NP <- nrow(V)
    weights <- rep(1, NP)
    if (prec_weights=='raw'){
        weights <- 1/diag(V)
        weights <- normalize_weights(x=weights)
    }
    return(weights)
}
