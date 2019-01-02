## File Name: linking_proc_itempars.R
## File Version: 0.07

linking_proc_itempars <- function(itempars)
{
    #*** convert itempars to data frame
    itempars <- as.data.frame( itempars )
    weights_exist <- TRUE
    # include wgt if there does not exist a fifth columm
    if ( ncol(itempars)==4){
        itempars$wgt <- 1
        weights_exist <- FALSE
    }
    # extract studies
    studies <- sort( paste( unique( itempars[,1] ) ) )
    NS <- length(studies)
    # extract items
    items <- sort( paste( unique( itempars[,2] ) ) )
    NI <- length(items)
    # define a and b matrices
    wgtM <- bM <- aM <- matrix(NA, nrow=NI, ncol=NS)
    rownames(wgtM) <- rownames(bM) <- rownames(aM) <- items
    colnames(wgtM) <- colnames(bM) <- colnames(aM) <- studies
    # define item parameters
    for (ss in studies){
        itempars.ss <- itempars[ itempars[,1]==ss, ]
        items_ss <- paste(itempars.ss[,2])
        aM[ items_ss, ss ] <- itempars.ss[,3]
        bM[ items_ss, ss ] <- itempars.ss[,4]
        wgtM[ items_ss, ss ] <- itempars.ss[,5]
    }
    wgtM <- wgtM / matrix( rowSums( wgtM, na.rm=TRUE ), nrow=NI, ncol=NS )
    est_pars <- ! is.na(wgtM)

    #--- output
    res <- list(itempars=itempars, NS=NS, NI=NI, items=items, studies=studies,
                wgtM=wgtM, aM=aM, bM=bM, est_pars=est_pars, weights_exist=weights_exist)
    return(res)
}
