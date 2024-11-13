## File Name: linking_haberman_lq_pw_le_grad.R
## File Version: 0.101


# define analytical gradient
linking_haberman_lq_pw_le_grad <- function(par_delta, par_gamma, des, h=1e-4)
{
    I <- des$I
    G <- des$G
    G1 <- G-1
    val <- matrix(0,nrow=I,ncol=2)
    design <- des$design
    w <- des$w
    ND <- nrow(design)
    coef_A <- c(1, par_delta[1L:(G-1)] )
    coef_B <- c(0, par_delta[(G-1)+1L:(G-1)] )
    NPD <- length(par_delta)
    grad <- matrix(0, nrow=I, ncol=NPD)

    for (dd in 1L:ND){
        gg_dd <- design$study1[dd]
        hh_dd <- design$study2[dd]
        ii_dd <- design$item[dd]

        #-- item discriminations
        w1 <- log( par_gamma[ 2*I*(gg_dd-1) + 2*(ii_dd-1)+1] )
        w2 <- log( par_gamma[ 2*I*(hh_dd-1) + 2*(ii_dd-1)+1] )
        e1 <- w[dd]*(w1-w2-log(coef_A[gg_dd])+log(coef_A[hh_dd]) )

        if (gg_dd>1){
            grad[ii_dd,gg_dd-1] <- grad[ii_dd,gg_dd-1] - e1
        }
        if (hh_dd>1){
            grad[ii_dd,hh_dd-1] <- grad[ii_dd,hh_dd-1] + e1
        }

        #-- item intercepts
        y1 <- par_gamma[ 2*I*(gg_dd-1) + 2*(ii_dd-1)+2] * coef_A[ gg_dd ]
        y2 <- par_gamma[ 2*I*(hh_dd-1) + 2*(ii_dd-1)+2] * coef_A[ hh_dd ]
        e2 <- w[dd]*(y1-y2-coef_B[gg_dd]+coef_B[hh_dd])
        if (gg_dd>1){
            grad[ii_dd,G1+gg_dd-1] <- grad[ii_dd,G1+gg_dd-1] - e2
        }
        if (hh_dd>1){
            grad[ii_dd,G1+hh_dd-1] <- grad[ii_dd,G1+hh_dd-1] + e2
        }
    }
    # grad <- grad/I

    return(grad)
}
