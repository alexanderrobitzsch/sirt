## File Name: rasch_mml2_calc_prob.R
## File Version: 0.06

rasch_mml2_calc_prob <- function( theta.k, b, fixed.a, fixed.c, fixed.d,
    alpha1, alpha2, Qmatrix, eps=1e-40)
{
    pjk <- prob_genlogis_4pl(theta=theta.k, b=b, a=fixed.a, c=fixed.c, d=fixed.d,
                        alpha1=alpha1, alpha2=alpha2, Qmatrix=Qmatrix)
    pjk <- ( pjk + eps ) / ( 1 + 2*eps )
    pjk.M <- t(pjk)
    return(pjk.M)
}
