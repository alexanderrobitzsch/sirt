## File Name: xxirt_nr_pml_preproc_data.R
## File Version: 0.11


xxirt_nr_pml_preproc_data <- function(em_args, pml_args)
{

    I <- em_args$I
    ncat <- em_args$ncat
    items <- em_args$items
    G <- em_args$G
    group <- em_args$group
    Theta <- em_args$Theta
    TP <- nrow(Theta)
    K <- max(ncat)
    if (is.null(pml_args$W1)){
        pml_args$W1 <- rep(1/I, I)
    }
    if (is.null(pml_args$W2)){
        W2 <- matrix( 2/I/(I-1), nrow=I, ncol=I)
        rownames(W2) <- colnames(W2) <- items
        W2[ upper.tri(W2)] <- 0
        diag(W2) <- 0
        pml_args$W2 <- W2
    }
    W2 <- pml_args$W2
    # convert W2 into a long format
    W2_long <- data.frame( item1=rep(1L:I,each=I), item2=rep(1L:I, I),
                                    w2=as.numeric(W2) )
    W2_long <- W2_long[ W2_long$w2 > 0, ]
    NI2 <- nrow(W2_long)
    dat1_resp <- em_args$dat1_resp
    weights <- em_args$weights

    freq1 <- array(0, dim=c(I,K,G))
    for (gg in 1L:G){
        for (ii in 1L:I){
            for (hh in 1L:K){
                freq1[ii,hh,gg] <- sum( dat1_resp[,ii,hh] * weights * (group==gg) )
            }
        }
    }
    freq2 <- array(0, dim=c(NI2, K, K,G) )
    nn <- 1
    for (gg in 1L:G){
        for (nn in 1L:NI2){
            for (hh in 1L:K){
                for (jj in 1L:K){
                    freq2[nn,hh,jj,gg] <- sum( dat1_resp[, W2_long[nn,1], hh ] *
                                    dat1_resp[, W2_long[nn,2], jj ] * weights *
                                            ( group==gg) )
                }
            }
        }
    }
    pml_args$W2_long <- W2_long
    pml_args$I <- I
    pml_args$G <- G
    pml_args$K <- K
    pml_args$TP <- TP
    pml_args$NI2 <- NI2
    pml_args$freq1 <- freq1
    pml_args$freq2 <- freq2

    #-- output
    return(pml_args)
}
