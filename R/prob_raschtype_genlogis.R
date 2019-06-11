## File Name: prob_raschtype_genlogis.R
## File Version: 0.11


# probability raschtype models
prob_raschtype_genlogis <- function( theta, b, alpha1, alpha2, fixed.a=1+0*b,
        Qmatrix=NULL )
{
    if ( is.null(Qmatrix) ){
        XX <- TAM::tam_outer(x=theta, y=b, op="-")
        XX <- as.vector(XX)
        aM <- sirt_matrix2(x=fixed.a, nrow=length(theta))
        XX <- aM * XX
    }
    if ( ! is.null(Qmatrix) ){
        XX0 <- tcrossprod( as.matrix(theta), Qmatrix )
        XX <- XX0 - outer( rep(1,nrow(theta)), b )
        XX <- as.vector(XX)
        aM <- sirt_matrix2(x=fixed.a, nrow=length(theta))
        XX <- aM * XX
    }
    pm <- pgenlogis(x=XX, alpha1=alpha1, alpha2=alpha2 )
    pm <- matrix( pm, ncol=length(b))
    return(pm)
}


.prob.raschtype.genlogis <- prob_raschtype_genlogis
