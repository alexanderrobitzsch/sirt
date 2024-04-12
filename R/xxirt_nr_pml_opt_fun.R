## File Name: xxirt_nr_pml_opt_fun.R
## File Version: 0.127

xxirt_nr_pml_opt_fun <- function(x, em_args, output_all=FALSE, use_rcpp=TRUE)
{

    Theta <- em_args$Theta
    pml_args <- em_args$pml_args

    #*** compute prior distribution
    prior_Theta <- xxirt_compute_prior_Theta_from_x(x=x, em_args=em_args)

    #* compute item response probabilities
    probs_items <- xxirt_compute_prob_item_from_x(x=x, em_args=em_args)
    eps <- 1e-14

    I <- pml_args$I
    K <- pml_args$K
    G <- pml_args$G
    TP <- pml_args$TP
    W1 <- pml_args$W1
    freq1 <- pml_args$freq1
    NI2 <- pml_args$NI2
    W2_long <- as.matrix(pml_args$W2_long)
    freq2 <- pml_args$freq2

    eval_fun_args <- list(prior_Theta=prior_Theta, probs_items=probs_items, freq1=freq1,
                            freq2=freq2, W1=W1, W2_long=W2_long, G=G, K=K,
                            I=I, TP=TP, NI2=NI2, eps=eps )

    if (use_rcpp){
        fun_pml_eval <- sirt_rcpp_xxirt_nr_pml_opt_fun
    } else {
        fun_pml_eval <- xxirt_nr_pml_opt_fun_R
    }

    res <- do.call( what=fun_pml_eval, args=eval_fun_args)
    res <- -res$val

    #-- output
    return(res)
}

