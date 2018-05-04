## File Name: rm_facets_est_b_rater.R
## File Version: 0.15


#########################################
# estimation of rater severity
rm_facets_est_b_rater <- function( b.item , b.rater , Qmatrix , tau.item ,
        VV , K , I , TP , a.item , a.rater , item.index , rater.index ,
        n.ik , numdiff.parm=.001 , max.b.increment=1 , theta.k , msteps ,
        mstepconv , b.rater.center , b.rater.fixed  )
{
    h <- numdiff.parm
    diffindex <- rater.index
    RR <- length(b.rater)    
    cat("  M steps b.rater parameter  |")
    it <- 0
    conv1 <- 1000
    #--- input for rm_facets_calcprobs
    args <- list( b.item=b.item, b.rater=b.rater, Qmatrix=Qmatrix, tau.item=tau.item, VV=VV, K=K, I=I, TP=TP, 
                    a.item=a.item, a.rater=a.rater, item.index=item.index, rater.index=rater.index, 
                    theta.k=theta.k, RR=RR )                         
    #---- begin algorithm
    while( ( it < msteps ) & ( conv1 > mstepconv ) ){
        b0 <- b.rater
        args$b.rater <- b.rater
        pjk <- do.call( what=rm_facets_calcprobs, args=args)
        args$b.rater <- b.rater + h
        pjk1 <- do.call( what=rm_facets_calcprobs, args=args)
        args$b.rater <- b.rater - h
        pjk2 <- do.call( what=rm_facets_calcprobs, args=args)
        #-- increment        
        res <- rm_numdiff_index( pjk=pjk, pjk1=pjk1, pjk2=pjk2, n.ik=n.ik, diffindex=diffindex, 
                    max.increment=max.b.increment, numdiff.parm=numdiff.parm ) 
        increment <- res$increment
        b.rater <- b.rater + increment
        if ( ! is.null( b.rater.fixed) ){        
            ind <- which( ! is.na( b.rater.fixed  ) )
            b.rater[ind] <- b.rater.fixed[ind]
            res$d2[ ind ] <- -1E10            
        }
        brc <- mean( b.rater )

        #-- centering
        b.rater <- rm_center_vector( vec=b.rater , center_type=b.rater.center)    
        conv1 <- max( abs( b.rater - b0 ) )
        it <- it+1
        cat("-")  
    }
    cat(" " , it , "Step(s) \n")    
    res <- list(b.rater = b.rater , se.b.rater = sqrt( abs(-1/res$d2 ) ) , 
                            ll = sum(res$ll0), brc = brc     )
    return(res)
}    


.rm.facets.est.b.rater <- rm_facets_est_b_rater
