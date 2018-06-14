## File Name: prob_genlogis_4pl.R
## File Version: 0.02

prob_genlogis_4pl <- function(theta, b, a, c, d, alpha1, alpha2, Qmatrix)
{
    pjk <- .prob.raschtype.genlogis( theta=theta, b=b, alpha1=alpha1,
                alpha2=alpha2, fixed.a=a, Qmatrix=Qmatrix)
    np <- nrow(pjk)
    cM <- sirt_matrix2( x=c, nrow=np )
    dM <- sirt_matrix2( x=d, nrow=np )
    pjk <- cM + (dM - cM) * pjk
    #--- OUTPUT
    return(pjk)
}
