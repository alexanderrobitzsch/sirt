## File Name: dif.strata.variance.R
## File Version: 0.151



#* stratified DIF variance
dif.strata.variance <- function( dif, se.dif, itemcluster )
{
    # stratified DIF variance
    # means in differential item functioning
    # itemcluster is a vector of strata corresponding to items
    stratadif <- stats::aggregate( 1+0*dif, list(itemcluster), sum, na.rm=TRUE )
    colnames(stratadif) <- c('strata', 'N.Items' )
    stratadif <- data.frame(stratadif)
    stratadif$M.DIF <- stats::aggregate( dif, list(itemcluster), mean, na.rm=TRUE )[,2]
    # DIF in strata
    SS <- nrow(stratadif)
    for (ss in 1:SS){
        items.ss <- which( itemcluster==stratadif[ss,'strata']  )
        dif.ss <- dif[ items.ss ]
        difv.ss <- dif.variance( dif=dif.ss, se.dif=se.dif[ items.ss ] )
        stratadif$weighted.tau[ss] <- difv.ss$weighted.DIFSD
        stratadif$unweighted.tau[ss] <- difv.ss$unweighted.DIFSD
    }
    stratadif[ is.na(stratadif ) ] <- 0

    sd_ni1 <- stratadif$N.Items-1
    weighted.DIFSD <- sum(stratadif$N.Items/sum(stratadif$N.Items)*stratadif$weighted.tau)
    unweighted.DIFSD <- sum( sd_ni1/ (sum(stratadif$N.Items)-1)*stratadif$unweighted.tau)

    #-- output
    res <- list( stratadif=stratadif, weighted.DIFSD=weighted.DIFSD,
                    unweighted.DIFSD=unweighted.DIFSD)
    return(res)
}


