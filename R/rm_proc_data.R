## File Name: rm_proc_data.R
## File Version: 0.57

##########################################
# Data preprocessing rater models
rm_proc_data <- function( dat, pid , rater, rater_item_int=FALSE, reference_rater=NULL )
{
    #--- define reference rater
    if ( is.null(reference_rater) ){
         reference_rater0 <- sort( paste( rater ))[1]
    }

    #--- item-specific rater parameters and rearrange dataset
    if (rater_item_int){
        if (is.null(reference_rater)){
            reference_rater <- reference_rater0
        }
        res <- rm_proc_create_pseudoraters( dat=dat, rater=rater, pid=pid,
                    reference_rater=reference_rater )
        dat <- res$dat
        rater <- res$rater
        pid <- res$pid
        reference_rater <- res$reference_rater
    }

    #-- create rater indices
    rater <- paste(rater)
    # create table of rater indizes
    rater.index <- data.frame( "rater" = sort( unique( rater )) )
    rater.index$rater.id <- seq( 1 , nrow(rater.index) )
    RR <- nrow(rater.index)

    # create table of person indizes
    person.index <- data.frame( "pid" = sort( unique( pid )) )
    person.index$person.id <- seq( 1 , nrow(person.index) )
    PP <- nrow( person.index )

    # number of variables
    VV <- ncol(dat)
    vars <- colnames(dat)
    # create data frame with crossed items and raters
    dat2 <- data.frame( matrix( NA , nrow=PP , ncol=RR*VV ) )
    colnames_dat2 <- colnames(dat2) <- paste0( rep(vars , RR ) , "-" ,
                rep( rater.index$rater , each=VV) )
    rownames(dat2) <- person.index$pid

    # create expanded dataset
    rater0 <- match( rater, rater.index$rater) - 1
    pid0 <- match( pid, person.index$pid) - 1
    dat2 <- sirt_rcpp_rm_proc_expand_dataset(dat=as.matrix(dat), rater0=rater0,
                    pid0=pid0, N=PP, R=RR)
    colnames(dat2) <- colnames_dat2
    dat20 <- dat2

    # variable list
    dataproc.vars <- list( item.index = rep( 1:VV , RR ),
                            rater.index = rep(1:RR , each=VV ) )

    # arrange response data
    N <- nrow(dat2)
    K <- max(dat2, na.rm=TRUE)
    I <- ncol(dat2)
    res <- sirt_rcpp_rm_proc_datasets_na_indicators(dat=as.matrix(dat2), K=K )
    dat2 <- res$dat2
    dat2.resp <- res$dat_resp
    colnames(dat2) <- colnames_dat2
    dat2.ind.resp <- array( res$dat2_ind_resp , dim=c(N,I,K+1) )

    #--- output
    res <- list( dat2=dat2, dat2.resp=dat2.resp, dat2.NA=dat20, dat=dat,
                person.index=person.index, rater.index=rater.index, VV=VV, N=PP, RR=RR,
                dataproc.vars=dataproc.vars, dat2.ind.resp=dat2.ind.resp, rater=rater, pid=pid, dat=dat,
                reference_rater=reference_rater)
    return(res)
}
#################################################

rm_proc <- rm_proc_data
.prep.data.rm <- rm_proc
