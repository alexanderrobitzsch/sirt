## File Name: xxirt_nr_opt_fun_pml_casewise.R
## File Version: 0.03


#-- case-wise optimization function for PML
xxirt_nr_opt_fun_pml_casewise <- function(x, object, eps=1e-14)
{
    em_args <- object$em_args
    pml_args <- em_args$pml_args
    prior_Theta <- xxirt_compute_prior_Theta_from_x(x=x, em_args=em_args)
    probs_items <- xxirt_compute_prob_item_from_x(x=x, em_args=em_args)
    I <- pml_args$I
    K <- pml_args$K
    G <- pml_args$G
    TP <- pml_args$TP
    W1 <- pml_args$W1
    freq1 <- pml_args$freq1
    NI2 <- pml_args$NI2
    W2_long <- as.matrix(pml_args$W2_long)
    freq2 <- pml_args$freq2
    group0 <- object$group - 1
    weights <- object$weights
    dat1 <- object$dat1
    dat_resp <- object$dat_resp==1
    eval_fun_args <- list(prior_Theta=prior_Theta, probs_items=probs_items,
                            W1=W1, W2_long=W2_long, G=G, K=K,
                            I=I, TP=TP, NI2=NI2, eps=eps, group0=group0,
                            weights=weights, dat1=dat1, dat_resp=dat_resp)
    res <- do.call( what=sirt_rcpp_xxirt_nr_pml_casewise_opt_fun, args=eval_fun_args)
    return(-res$case_val)
}
