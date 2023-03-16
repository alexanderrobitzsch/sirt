## File Name: gom_em_item_parameters.R
## File Version: 0.09


gom_em_item_parameters <- function(dat2, dat2.resp, model, b, lambda, K,
    weights, progress)
{
    item <- data.frame("item"=colnames(dat2))
    item$N <- colSums( weights*dat2.resp )
    item$p <- colSums( weights*dat2, na.rm=TRUE) / item$N
    item$b <- b
    if (model !="GOMRaschxxx"){
        for (kk in 1:K){
            item[,paste0("lam.Cl",kk)] <- lambda[,kk]
        }
    }
    obji <- item
    for (vv in seq(2,ncol(obji) )){
        obji[,vv] <- round( obji[,vv], 3 )
    }
    if (progress){
        cat("*********************************\n")
        cat("Item Parameters\n")
        print(obji)
    }
    #-- output
    return(item)
}
