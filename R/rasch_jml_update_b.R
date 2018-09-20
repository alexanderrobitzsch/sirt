## File Name: rasch_jml_update_b.R
## File Version: 0.02



# update item difficulty estimation (Rasch model)
rasch_jml_update_b <- function( b, theta, freq.thetapattern, freq.dat.resp.thetapattern,
            constraints=NULL, conv=.0001, suffB, progress=progress, bsteps=4)
{
    b.change <- 1
    iter <- 0
    while( max( abs( b.change  ) ) > conv & ( iter < bsteps ) ){
        p.ia <- stats::plogis( theta, matrix( b, nrow=length(theta), length(b), byrow=T ) )
        deriv <- colSums( - freq.thetapattern * p.ia * ( 1- p.ia ) )
        diff <- suffB + colSums( freq.dat.resp.thetapattern *  p.ia  )
        b.change <-  diff / deriv
        if (! is.null(constraints)){
            b.change[ constraints[,1] ] <- 0
        }
        if (progress){
            cat("-")
            utils::flush.console()
        }
        b <- b - b.change
        iter <- iter + 1
    }
    return(b)
}


.update.b.rasch.jml2 <- rasch_jml_update_b
