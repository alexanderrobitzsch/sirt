## File Name: summary_round_helper.R
## File Version: 0.08

summary_round_helper <- function( obji, digits, exclude=NULL, print=TRUE){
    NC <- ncol(obji)
    ind <- 1:NC
    if ( ! is.null(exclude) ){
        ind2 <- which( colnames(obji) %in% exclude )
        ind <- setdiff( ind, ind2 )
    }
    obji[,ind] <- round( obji[,ind], digits )
    rownames(obji) <- NULL
    print(obji)
    invisible(obji)
}
