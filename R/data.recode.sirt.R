## File Name: data.recode.sirt.R
## File Version: 0.02


#*** utility function for recoding a raw dataset
data.recode.sirt <- function( data.raw, keys )
{
    item.stat <- keys
    V <- ncol(data.raw)
    data.scored <- matrix( 0, nrow(data.raw), ncol(data.raw) )
    colnames(data.scored) <- colnames(data.raw )
    for (vv in 1L:V){
        data.scored[,vv] <- 1* ( paste(data.raw[,vv])==
                    paste(item.stat[ item.stat$item==colnames(data.raw)[vv], 'key' ]) )
        data.scored[ paste( data.raw[,vv] )=='NA', vv ] <- NA
    }
    return(data.scored)
}
