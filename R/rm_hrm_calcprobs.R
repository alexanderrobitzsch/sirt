## File Name: rm_hrm_calcprobs.R
## File Version: 0.08

################################################################
# calculate probabilities
rm_hrm_calcprobs <- function(  c.rater , Qmatrix , tau.item ,
                VV , K , I , TP , a.item , d.rater , item.index , rater.index ,
                theta.k ,RR , prob.item=NULL , prob.rater = NULL, output_prob_total=FALSE )
{
    # calculate probabilities for true ratings
    a <- a.item
    b <- tau.item
    if ( is.null( prob.item ) ){ 
        res <- .rm.pcm.calcprobs( a=a, b=b, Qmatrix=Qmatrix, theta.k=theta.k, I=VV, K=K, TP=TP ) 
    } else { 
        res <- prob.item 
    }
    # calculate probabilities for raters
    calc.rater <- FALSE
    if (is.null(prob.rater)){
        calc.rater <- TRUE
    }
    dimA <- c(I , K+1, K+1 )
    res2 <- res[ item.index ,,]    
    dimB <- dim(res2)
    BM <- matrix( res2 , dimA[1]*dimB[2] , dimB[3] )
    #****    
    # if prob.rater is calculated    
    if (calc.rater){
        res2 <- rm_sdt_probraterfct1( crater=c.rater, drater=d.rater, dimA=dimA, B=BM, dimB=dimB ) 
        prob.categ <- array( res2$probtotal , dim= c(dimA[c(1,2)],dimB[3]) )
        prob.rater <- array( res2$PRA , dim=dimA )        
    }
    #***
    # if prob.rater is not calculated
    AM <- matrix( prob.rater , dimA[1]*dimA[2] , dimA[3] )
    if ( ! calc.rater ){
        y <- rm_sdt_arraymult1( A=AM, dimA=dimA, B=BM, dimB=dimB )
        prob.categ <- array( y , dim= c(dimA[c(1,2)],dimB[3]) )    
    }
    #-- restrict prob.total
    prob.total <- prob.categ
    eps <- 1E-4
    prob.total[ prob.total < eps ] <- eps    
    #--- output
    res <- list(prob.total=prob.total , prob.rater=prob.rater, prob.item=res )
    if (output_prob_total){
        res <- res$prob.total
    }
    return(res)    
}
#############################################################

.rm.hrm.calcprobs <- rm_hrm_calcprobs
