## File Name: qmc.nodes.R
## File Version: 0.06
## File Last Change: 2017-01-18 11:02:51
qmc.nodes <- function( snodes , ndim ){          
    r1 <- sfsmisc::QUnif(n=snodes, min = 0, max = 1, n.min = 1, p=ndim, leap = 409)                     
    theta <- as.matrix( stats::qnorm( r1 ) )
    return(theta)
        }
