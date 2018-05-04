## File Name: rm_hrm_est_d_rater.R
## File Version: 0.14

###################################################            
# d.rater
rm_hrm_est_d_rater <- function(  c.rater , Qmatrix , tau.item ,
                    VV , K , I , TP , a.item , d.rater , item.index , rater.index ,
                    n.ik , numdiff.parm, max.b.increment=1,theta.k ,
                    msteps , mstepconv , d.min , d.max , est.d.rater , prob.item, d.rater0, diffindex,
                    d.prior )
{
    h <- numdiff.parm            
    
    RR <- I/VV    
    cat("  M steps d.rater parameter    |")
    it <- 0
    conv1 <- 1000
    ND <- length(d.rater)
    
    #-- input calcprobs
    args <- list( c.rater=c.rater, Qmatrix=Qmatrix, tau.item=tau.item, VV=VV, K=K, I=I, TP=TP, a.item=a.item, 
                    d.rater=d.rater, item.index=item.index, rater.index=rater.index, theta.k=theta.k, RR=RR, 
                    prob.item=prob.item, prob.rater=NULL, output_prob_total=TRUE )

    #--- begin M-steps
    while( ( it < msteps ) & ( conv1 > mstepconv ) ){
        b0 <- d.rater

        args$d.rater <- d.rater
        pjk <- do.call(what=rm_hrm_calcprobs, args=args)
        
        args$d.rater <- d.rater + h
        pjk1 <- do.call(what=rm_hrm_calcprobs, args=args)

        args$d.rater <- d.rater - h
        pjk2 <- do.call(what=rm_hrm_calcprobs, args=args)
        
        #-- increments
        res <- rm_numdiff_index( pjk=pjk, pjk1=pjk1, pjk2=pjk2, n.ik=n.ik, diffindex=diffindex, 
                    max.increment=max.b.increment, numdiff.parm=numdiff.parm, prior=d.prior, 
                    value=d.rater )     
        d.rater <- d.rater + res$increment[diffindex]
        d.rater[ d.rater < d.min ] <- d.min        
        d.rater[ d.rater > d.max ] <- d.max                
        #max.b.increment <- abs( b.rater - b0 )
        conv1 <- max( abs( d.rater - b0 ) )
        it <- it+1
        cat("-")  
    }
    
    #- trim increments
    d.rater <- rm_trim_increments_mstep( parm=d.rater, parm0=d.rater0 , max.increment=max.b.increment )    
    
    cat(" " , it , "Step(s) \n")    
    res <- list(d.rater = d.rater , se.d.rater = sqrt( abs(-1/res$d2) ) , ll = sum(res$ll0) )
    return(res)
}                

.rm.hrm.est.d.rater <- rm_hrm_est_d_rater
            
            
