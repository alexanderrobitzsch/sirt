## File Name: noharm_sirt_preproc_pattern_matrix.R
## File Version: 0.12

noharm_sirt_preproc_pattern_matrix <- function(pattmat, minval=0, symm=FALSE)
{
    I <- nrow(pattmat)
    D <- ncol(pattmat)
    if (symm){
        pattmat <- ( pattmat + t(pattmat) ) / 2
    }
    mp <- max(pattmat)
    ind <- ( pattmat > 0 ) & ( pattmat <=1 )
    if (sum(ind) > 0 ){
        DI <- prod(dim(pattmat))
        v1 <- matrix( 1:DI, nrow=I )
        if (symm){
            v1[upper.tri(v1)] <- v1[lower.tri(v1)]
        }
        pattmat[ind] <- mp + v1[ind]
    }
    pattmat[ pattmat==0] <- NA
    elem <- stats::na.omit(as.vector(pattmat))
    if (length(elem)>0){
        pattmat <- matrix( match( pattmat, unique(elem) ), nrow=I) + minval
    }
    return(pattmat)
}

