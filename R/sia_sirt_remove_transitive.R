## File Name: sia_sirt_remove_transitive.R
## File Version: 0.13
## File Last Change: 2018-12-30


#**** remove transitive relations
sia_sirt_remove_transitive <- function(I1){
    I <- ncol(I1)
    diag(I1) <- 0
    BB <- 1
    IS <- I1
    iter <- 0
    while (BB > .0001){
        I0 <- IS
        iter <- iter + 1
        for (ii in 1:I ){
            for (jj in 1:I){
                if (ii!=jj ){
                    if ( sum( IS[ii,] * IS[,jj] ) > 0 ){
                            I1[ii,jj] <- 0
                    }
                }
            }
        }
        IS <- IS + IS %*% I1
        IS <- 1*(IS > 0 )
        BB <- sum( abs( I0 - IS ) )
    }
    return(I1)
}


.sia.remove.transitive <- sia_sirt_remove_transitive
