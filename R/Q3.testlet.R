## File Name: Q3.testlet.R
## File Version: 1.21
## File Last Change: 2018-12-30


################################################################################
# Summarizing testlet effect using Q3 statistic
Q3.testlet <- function( q3.res, testlet.matrix, progress=TRUE )
{
    N.item.testlet <- stats::aggregate( rep(1, nrow(testlet.matrix) ),
                              list( testlet.matrix[,1]), sum )
    testlet.matrix <- testlet.matrix[ testlet.matrix[,1] %in%
                                N.item.testlet[ N.item.testlet[,2] > 1, 1 ], ]
    testlets <- sort( unique( testlet.matrix[,1] ) )
    testlet.q3 <- t( sapply( testlets, FUN=function(testlet){
            testlet.items <- testlet.matrix[ testlet.matrix[,1]==testlet, 2 ]
            ti.ind <- colnames(q3.res$q3.matrix) %in% testlet.items
            # c( sum( ti.ind), mean( q3.res$q3.matrix[ ti.ind, ti.ind ], na.rm=T ) )
            # correction thanks to Thomas Kiefer (2014-03-06)
                c( sum( ti.ind), mean(
                     q3.res$q3.matrix[ ti.ind, ti.ind ][
                            lower.tri( diag( sum(ti.ind))) ],
                     na.rm=TRUE ) )
                } ) )
    colnames(testlet.q3) <- c("N.Items", "Mean.Q3" )
    testlet.q3 <- data.frame( "Testlet"=testlets, testlet.q3,
                            "mean"=mean(q3.res$q3.long[,3]) )
    rownames(testlet.q3) <- NULL
    # mean Q3-statistics between testlets
    TT <- length(testlets)
    matr <- matrix( 1, nrow=TT, ncol=TT )
    colnames(matr) <- rownames(matr) <- testlets
    for (ii1 in seq(1,TT-1)){
        for (ii2 in seq(ii1+1, TT )){
            tt1 <- paste(testlets[ii1])
            tt2 <- paste(testlets[ii2])
            itt1 <- testlet.matrix[ testlet.matrix[,1]==tt1,2 ]
            itt2 <- testlet.matrix[ testlet.matrix[,1]==tt2,2 ]
            q.tt <- q3.res$q3.matrix[ paste( itt1 ), paste( itt2) ]
            matr[ tt1, tt2 ] <- matr[tt2,tt1] <- mean( q.tt, na.rm=TRUE )
        }
    }
    diag(matr) <- testlet.q3$Mean.Q3
    if (progress){
        cat( "\nMean Q3 Testlets:", round( mean(testlet.q3$Mean.Q3), 5 ),"\n\n")
        print( testlet.q3, digits=3 )
        cat( "\n\nMean Q3 between testlets \n\n")
        matr1 <- round( matr, 3 )
        print( matr1, digits=3 )
    }
    #--- OUTPUT
    res <- list( "testlet.q3"=testlet.q3, "testlet.q3.korr"=matr )
    return(res)
}
##########################################################################

# Q3.testlet <- testlet.yen.q3
