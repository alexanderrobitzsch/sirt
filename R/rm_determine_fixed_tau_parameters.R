## File Name: rm_determine_fixed_tau_parameters.R
## File Version: 0.05

rm_determine_fixed_tau_parameters <- function( K, maxK, VV, tau.item.fixed=NULL, val=99)
{
    if ( min(maxK) < K ){
        for (vv in 1:VV){
            K.vv <- maxK[vv]
            if ( K.vv < K ){
                for (zz in (K.vv+1):K ){
                    d1 <- data.frame( "item"=vv, "categ"=zz, "val"=val)
                    tau.item.fixed <- rbind( tau.item.fixed, d1 )
                }
            }
        }
        tau.item.fixed <- as.matrix(tau.item.fixed )
    }
    return(tau.item.fixed)
}
