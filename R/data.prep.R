## File Name: data.prep.R
## File Version: 1.153

#----- data preparations for rasch.jml and rasch.mml
data.prep <- function( dat, weights=NULL, use.freqpatt=TRUE,
        standardize_weights=TRUE)
{
    item.means <- colMeans( dat, na.rm=TRUE )
    item.elim <- which( item.means %in% c(0,1))
    if ( length( item.elim ) > 0 ){
        stop( cat( paste( 'There are', length(item.elim), 'Items with no variance!') ) )
    }
    if ( any( is.na(item.means)) ){
        stop( 'There are items which contains only missings!')
    }
    n <- nrow(dat)
    I <- ncol(dat)
    if( is.null(weights) ){  weights <- rep( 1, n ) }
    # indicator for nonmissing response
    dat.9 <- dat
    dat.9[ is.na(dat) ] <- 9

    #* pattern
    if ( use.freqpatt ){
        freq.patt <- apply( dat.9, 1, FUN=function(ll){ paste(ll, collapse='' ) } )  #
        dat1 <- data.frame( table( freq.patt ) )
    } else {
        freq.patt <- paste('FP', 1000000 + 1L:n, sep='')
        dat1 <- data.frame( freq.patt )
        colnames(dat1)[1] <- 'freq.patt'
    }
    # weighting the frequencies if survey weights are supplied
    if( !is.null(weights) ){
        # standardize weights
        if (standardize_weights){
            weights <- weights / sum(weights) * n
        }
        if ( use.freqpatt ){
            dat1[,2] <- stats::aggregate( weights, list( freq.patt), sum )[,2]
        } else {
            dat1[,'Freq'] <- weights
        }
    }
    # item pattern corresponding to frequency pattern
    if ( use.freqpatt){
        dat2 <- matrix( as.numeric( unlist( strsplit( paste(dat1[,1]), '' ) ) ),
                            ncol=ncol(dat), byrow=TRUE)
    } else {
        dat2 <- dat.9 }
        dat2.resp <- 1*( dat2 !=9 )
        dat2[ dat2==9 ] <- 0
        # mean right
        dat1$mean <- rowSums( dat2 * dat2.resp )  / rowSums( dat2.resp )
        freq.patt <- data.frame(  freq.patt, rowMeans( dat, na.rm=TRUE ), 1L:n )
        colnames(freq.patt)[2L:3] <- c('mean', 'index' )
        list( dat=dat, dat2=dat2, dat2.resp=dat2.resp, dat1=dat1,
                freq.patt=freq.patt, I=I, n=n, dat9=dat.9 )
}


.data.prep <- data.prep



#*** Small function which helps for printing purposes
.prnum <- function( matr, digits )
{
    VV <- ncol(matr)
    for (vv in 1L:VV){
        if ( is.numeric( matr[,vv]) ){
            matr[,vv] <- round( matr[,vv], digits )
        }
    }
    print(matr)
}


#-- Function for calculation of a response pattern
#   for dichotomous responses
resp.pattern2 <- function(x)
{
    n <- nrow(x)
    p <- ncol(x)
    mdp <- (x %*% (2^((1L:ncol(x)) - 1))) + 1
    misspattern <- mdp[,1]
    misspattern <- list( miss.pattern=mdp[,1],
                        mp.index=match( mdp[,1], sort( unique(mdp[,1] ) ) ) )
    return(misspattern)
}

