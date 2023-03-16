## File Name: rasch_jml_centerpersons.R
## File Version: 0.04
## File Last Change: 2018-12-30

rasch_jml_centerpersons <- function(theta, dat1, centerpersons)
{
    if (centerpersons){
        theta <- theta - stats::weighted.mean( theta, dat1[,2] )
    }
    return(theta)
}
