## File Name: monoreg.rowwise.R
## File Version: 0.06

##############################################################
# monotone regression for all rows in a matrix
monoreg.rowwise <- function(yM,wM){
    yM <- as.matrix(yM)
    wM <- as.matrix(wM)
    res <- monoreg_rowwise_Cpp( yM, wM )
    return(res)
                    }
##############################################################
# monotone regression for all columns in a matrix
monoreg.colwise <- function(yM,wM){
    yM <- as.matrix(yM)
    wM <- as.matrix(wM)
    res <- monoreg_rowwise_Cpp( t(yM), t(wM) )
    return(t(res))
}
