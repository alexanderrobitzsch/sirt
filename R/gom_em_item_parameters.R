## File Name: gom_em_item_parameters.R
## File Version: 0.02


gom_em_item_parameters <- function(dat2, dat2.resp, model, b, lambda, K, progress)
{
    item <- data.frame("item"=colnames(dat2))
    item$N <- colSums( dat2.resp )
    item$p <- colMeans( dat2, na.rm=TRUE)
    item$b <- b
    if (model !="GOMRaschxxx"){
        for (kk in 1:K){
            item[,paste0("lam.Cl",kk)] <- lambda[,kk]
        }
    }
    obji <- item
    for (vv in seq(2,ncol(obji) )){
        obji[,vv] <- round( obji[,vv],3 )
    }
    if (progress){
        cat("*********************************\n")
        cat("Item Parameters\n")
        print(obji)
    }
    #-- output
    return(item)
}
