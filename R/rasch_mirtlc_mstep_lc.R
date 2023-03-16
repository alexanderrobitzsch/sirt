## File Name: rasch_mirtlc_mstep_lc.R
## File Version: 0.04
## File Last Change: 2019-09-14


#---- calculate class probabilities
rasch_mirtlc_mstep_lc <- function( pjk, n.k, r.jk, n.jk, G, Nclasses )
{
    if (G==1){
        pi.k <- n.k / sum( n.k )
        pi.k <- matrix( pi.k, nrow=ncol(pjk), ncol=nrow(pjk) )
    }
    if ( G> 1){
        pi.k <- n.k / matrix( colSums(n.k ), nrow=Nclasses, ncol=G, byrow=TRUE)
    }
    for (cc in 1:Nclasses ){
        if (G==1){
            pjk[ cc, ] <- r.jk[, cc, 1] / n.jk[, cc, 1]
        }
        if (G>1){
            pjk[ cc, ] <- rowSums( r.jk[, cc, ] ) / rowSums( n.jk[, cc, ]  )
        }
    }
    res <- list( pi.k=pi.k, pjk=pjk )
    return(res)
}


.m.step.mirtlc.lc <- rasch_mirtlc_mstep_lc
