## File Name: write.fwf2.R
## File Version: 1.10


##############################################################################
write.fwf2 <- function( dat, format.full, format.round, savename ){
        if (is.null( colnames(dat) ) ){
            colnames(dat) <- paste( "V", 1:( ncol(dat) ), sep="")
                            }
        matr <- matrix( " ", nrow=nrow(dat), ncol=ncol(dat) )
        ind1 <- which( format.round <=0  )
        format.full[ ind1 ] <- format.full[ind1]
        format.round[ ind1 ] <- format.round[ind1]
        for (vv in 1:( ncol(matr) ) ){
            fvv <- format.round[vv]
            fff <- format.full[vv]
            matr[,vv] <- .write.format2( vec1=dat[,vv], ff=fff, fr=fvv )
                }
        matr <- apply( matr, 1, FUN=function(ll){ paste( ll,
                        collapse="" ) } )
        if ( is.vector(matr) ){
            writeLines( matr, paste( savename, ".dat", sep="") )
                    } else {
            utils::write.table( matr, paste( savename, ".dat", sep=""),
                        row.names=F, col.names=F)
                        }
        dfr <- data.frame( "variable"=colnames(dat),
                    "begin"=c( 1, cumsum( format.full )[ - ncol(dat) ]
                                + 1 ),
                    "end"=cumsum( format.full ),
                    "length"=format.full
                            )
        utils::write.table( dfr, paste( savename, "__LEGEND.txt",sep=""),
                quote=FALSE, row.names=FALSE, col.names=TRUE)

        return(dfr)
        }
##############################################################################

