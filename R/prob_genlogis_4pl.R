## File Name: prob_genlogis_4pl.R
## File Version: 0.08
## File Last Change: 2019-05-10


prob_genlogis_4pl <- function(theta, b, a, c, d, alpha1, alpha2, Qmatrix)
{
    pjk <- prob_raschtype_genlogis( theta=theta, b=b, alpha1=alpha1,
                alpha2=alpha2, fixed.a=a, Qmatrix=Qmatrix)
    if ( (any(c>0)) | (any(d<1)) ){
        np <- nrow(pjk)
        cM <- sirt_matrix2( x=c, nrow=np )
        dM <- sirt_matrix2( x=d, nrow=np )
        pjk <- cM + (dM - cM) * pjk
    }
    #--- output
    return(pjk)
}
