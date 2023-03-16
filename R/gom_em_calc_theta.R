## File Name: gom_em_calc_theta.R
## File Version: 0.17
## File Last Change: 2019-05-12


#--- calculate theta grid
gom_em_calc_theta <- function( K, problevels, eps=1e-5 )
{
    m1 <- problevels
    if ( ! is.matrix(problevels) ){
        PL <- length(problevels)
        m1 <- matrix(problevels, PL, 1 )
        for (kk in 2:K){
            NM <- nrow(m1)
            m1 <- cbind( m1[ rep( 1:NM, PL), ], rep( problevels, each=NM)  )
            m1 <- m1[ rowSums(m1) <=1, ]
        }
    }
    m1 <- m1[ abs( rowSums(m1) - 1 ) < eps, ]
    return(m1)
}



.gom.calc.theta <- gom_em_calc_theta
