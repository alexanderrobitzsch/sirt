## File Name: split_conc.R
## File Version: 0.01


#-- splits a string syn.ll and concatanates it with string vv
split_conc <- function( syn.ll, vv )
{
    g1 <- strsplit( syn.ll, vv, perl=FALSE )[[1]]
    Lg1 <- length(g1)
    vec <- NULL
    if (Lg1==1 ){
        vec <- c( g1, vv )
    }
    if (Lg1 > 1 ){
        vec <- rep("", Lg1 + (Lg1-1) )
        vec[ seq( 1, 2*Lg1, 2 ) ] <- g1
        vec[ seq( 2, 2*Lg1 - 1, 2 ) ] <- vv
        Ls1 <- nchar(syn.ll)
        if ( substring( syn.ll, Ls1, Ls1 )==vv ){
            vec <- c( vec, vv )
        }
    }
    return(vec)
}
