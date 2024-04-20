## File Name: linking_haberman_itempars_prepare.R
## File Version: 0.132


linking_haberman_itempars_prepare <- function(b, a=NULL, wgt=NULL)
{
    I <- nrow(b)
    NS <- ncol(b)
    if ( is.null(rownames(b) )){
        i0 <- ceiling( log10(I) + 1 )
        rownames(b) <- paste0('I', 10^i0 + 1L:I)
    }
    if ( is.null(colnames(b) )){
        colnames(b) <- 1L:NS
    }
    if (is.null(a)){
        a <- matrix(1, nrow=I, ncol=NS)
    }
    if (is.null(wgt)){
        wgt <- matrix(1, nrow=I, ncol=NS)
    }
    #-- arrange table
    itempars <- data.frame(study=rep(colnames(b), each=I),
                        item=rep(rownames(b), NS), a=as.vector(a),
                        b=as.vector(b), wgt=as.vector(wgt) )
    return(itempars)
}
