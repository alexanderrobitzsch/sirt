## File Name: stratified.cronbach.alpha.R
## File Version: 0.257


# stratified Cronbach's Alpha
stratified.cronbach.alpha <- function( data, itemstrata=NULL )
{
    stratcomp <- TRUE
    if ( is.null(itemstrata) ){
        itemstrata <- cbind( colnames(data), 1 )
        stratcomp <- FALSE
    }
    #---------------------------------------
    dfr <- data.frame( scale="total", stratified_cronbach_alpha_compute_alpha(data) )
    # calculation of stratified alpha
    itemstrata.u <- sort(unique( itemstrata[,2] ))
    for (gg in itemstrata.u){
        data_gg <- data[, itemstrata[ itemstrata[,2]==gg, 1] ]
        dfr1 <- data.frame( scale=gg,
            stratified_cronbach_alpha_compute_alpha(data=data_gg) )
        dfr <- rbind( dfr, dfr1 )
    }
    # stratified alpha
    dfr$alpha.stratified <- NA
    var_tot <- dfr[ -1, "var.tot" ]
    dfr_alpha <- dfr[ -1, "alpha" ]
    # dfr$alpha.stratified[1] <- 1 - sum (( 1 - dfr_alpha )*var_tot ) / var_tot
    dfr$alpha.stratified[1] <- 1 - sum (( 1 - dfr_alpha )*var_tot ) / dfr[ 1, "var.tot"]
    obji <- dfr
    obji[, -c(1:2)] <- round( obji[,-c(1:2) ], 3 )
    if ( ! stratcomp ){
        obji <- obji[1,]
    }
    print(obji)
    invisible(dfr)
}
