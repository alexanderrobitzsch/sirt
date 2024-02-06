## File Name: xxirt_postproc_parameters.R
## File Version: 0.236



xxirt_postproc_parameters <- function( partable, customTheta,
        items, probs_items, np_fun_item=NULL )
{
    #**** item parameters
    p1 <- partable[ partable$parfree==1, ]
    par_items <- p1$value
    names(par_items) <- p1$parlabel
    #*** theta distribution parameters
    cs <- customTheta
    par_Theta <- cs$par[ cs$est ]
    #*** structured form of parameters
    I <- length(items)
    parnames <- sort( unique( paste( partable$parname) ) )
    PN <- length(parnames)
    m1 <- matrix(NA, nrow=I, ncol=PN)
    rownames(m1) <- items
    colnames(m1) <- parnames
    for (pp in 1:PN){
        p1 <- partable[ partable$parname==parnames[pp], ]
        m1[ p1$item, parnames[pp] ] <- p1$value
    }
    p1 <- partable[ ! duplicated(partable$item ), ]
    dfr <- data.frame( item=items, type=paste(p1$type), m1  )
    rownames(dfr) <- NULL

    #*** probabilities item parameters
    pi_dim <- dim(probs_items)
    dimnames(probs_items)[[1]] <- items
    dimnames(probs_items)[[2]] <- paste0('Cat', seq(0,pi_dim[2]-1) )
    #*** lower and upper bounds
    p1 <- partable[ partable$parfree==1, c('item', 'type', 'parname',
                        'value', 'lower', 'upper', 'parindex', 'parlabel' ) ]
    p1$active <- 1 * ( p1$value > p1$lower )
    p1$active <- p1$active * ( p1$value < p1$upper )
    par_items_bounds <- p1

    np_item <- NULL
    if ( ! is.null(np_fun_item) ){
        np_item <- np_fun_item(x=par_items)
    }

    #*** output
    res <- list( par_items=par_items, par_Theta=par_Theta,
                    probs_items=probs_items, par_items_summary=dfr,
                    par_items_bounds=par_items_bounds, np_item=np_item )
    return(res)
}

