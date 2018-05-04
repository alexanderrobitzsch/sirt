## File Name: xxirt_proc_ParTable.R
## File Version: 0.36

#################################################
# process parameter table
xxirt_proc_ParTable <- function( itemtype , partable , items )
{
    #*** extract item types from partable
    itemtype <- unlist( sapply( items , FUN = function(ii){ 
                        partable$type[ paste(partable$item) == ii ][1] 
                        } )  )
    I <- length(items)
    partable$rowindex <- seq( 1 , nrow(partable) )
    #*** parameter index
    partable[ ! partable$est , "parindex" ] <- 0
    # indices <- sort( unique( partable$parindex ) )
    indices <- unique( partable$parindex )
    indices <- c( 0 , setdiff( indices , c(0) ) )
    IN <- length(indices)
    partable$parindex <- match( partable$parindex , indices ) - 1
    partable$parindex[ partable$parindex == 0 ] <- NA
    #*** set prior distributions of fixed parameters to NA
    partable[ ! partable$est , "prior" ] <- NA
    #*** list with parameter table indices        
    partable$parfree <- 1*partable$est
    partable_index <- as.list( 1:I )        
    for (ii in 1:I){
        partable_index[[ii]] <- which( partable$itemnr == ii )
    }
    ind <- which( duplicated( partable$parindex ) & ( ! is.na( partable$prior ) ) )    
    if ( length(ind) > 0 ){
        partable[ind , c("prior","prior_par1","prior_par2", "parlabel") ] <- NA        
    }    
    ind <- which( duplicated( partable$parindex ) )    
    if ( length(ind) > 0 ){
        partable[ind,"parfree"] <- 0
    }
    #*** extract ncat and maxK
    p1 <- partable[ ! duplicated( partable$item) , ]
    ncat <- p1$ncat    
    names(ncat) <- paste(p1$item)
    maxK <- max(ncat)
    #*** extract M-step method
    p1 <- partable[ partable$parfree == 1 , ]
    m1 <- TRUE
    if ( nrow(p1) > 0 ){
        m1 <- ( mean( p1$lower == - Inf ) < 1 ) | ( mean( p1$upper == Inf ) < 1 )
    }
    mstep_method <- if (m1){ "L-BFGS-B" } else { "BFGS" }
        
    #**** item indices per parameter
    NP <- -Inf
    if ( sum(partable$est) > 0 ){    
        NP <- max( partable$parindex , na.rm=TRUE )
    }        
    if ( NP > -Inf){
        item_index <- as.list( 1:NP )
        for (pp in 1:NP){
            p1 <- partable[ paste0( partable$parindex) == pp , ]
            names(item_index)[pp] <- p1[1,"parlabel"]
            item_index[[pp]] <- p1$itemnr
        }
    } else {
        item_index <- list()
    }
                
    #----- output
    res <- list( itemtype = itemtype , partable = partable ,
                    partable_index = partable_index , ncat = ncat ,
                    maxK = maxK , mstep_method = mstep_method , 
                    item_index=item_index )
    return(res)
}
################################################
