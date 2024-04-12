## File Name: xxirt_nr_pml_grad_fun.R
## File Version: 0.129

xxirt_nr_pml_grad_fun <- function(x, em_args, output_all=FALSE, use_rcpp=TRUE)
{

    Theta <- em_args$Theta
    pml_args <- em_args$pml_args
    h <- em_args$h

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

    NP <- em_args$NP
    NPI <- em_args$NPI
    NPT <- em_args$NPT

    #- define more arguments
    parindex_items <- em_args$parindex_items
    parindex_Theta <- em_args$parindex_Theta
    pml_args$in_par_Theta <- em_args$pml_args$in_par_Theta

    #- define gradient
    eval_fun_args <- list(prior_Theta=prior_Theta, probs_items=probs_items, freq1=freq1,
                            freq2=freq2, W1=W1, W2_long=W2_long, G=G, K=K,
                            I=I, TP=TP, NI2=NI2, eps=eps )
    fun_pml_eval <- sirt_rcpp_xxirt_nr_pml_opt_fun
    pml_fun0 <- do.call( what=fun_pml_eval, args=eval_fun_args)

    grad <- 0*x
    eval_grad_args <- eval_fun_args
    eval_grad_args$NP <- NP
    eval_grad_args$val1 <- pml_fun0$val1
    eval_grad_args$val2 <- pml_fun0$val2

    eval_grad_args$der_prior_Theta <- 0*prior_Theta
    eval_grad_args$der_probs_items <- der_probs_items0 <- 0*probs_items

    for (tpp in 1L:NP){
        pp_Theta <- pml_args$in_par_Theta[tpp]

        x1 <- sirt_add_increment(x=x, pos=tpp, value=h)
        if (pp_Theta){
            prior_Theta1 <- xxirt_compute_prior_Theta_from_x(x=x1, em_args=em_args)
            der_prior_Theta <- ( (prior_Theta1-prior_Theta)/h )
            eval_grad_args$der_prior_Theta <- der_prior_Theta
        } else {
            itpp <- em_args$item_index[[tpp]]
            probs_items1 <- xxirt_compute_prob_item_from_x(x=x1, em_args=em_args,
                                    item_index=itpp)
            der_probs_items <- der_probs_items0
            der_probs_items[itpp,,] <- ( (probs_items1-probs_items[itpp,,,drop=FALSE])/h )
            eval_grad_args$der_probs_items <- der_probs_items
        }

        eval_grad_args$pp_Theta <- pp_Theta
        eval_grad_args$index_freq1 <- em_args$pml_args$index_freq1[[tpp]]
        eval_grad_args$index_freq2 <- em_args$pml_args$index_freq2[[tpp]]

        grad[tpp] <- do.call( what=sirt_rcpp_xxirt_nr_pml_grad_fun_eval,
                                    args=eval_grad_args)

    }

    #-- output
    return(grad)
}

