## File Name: rasch_jml_centerpersons.R
## File Version: 0.02

rasch_jml_centerpersons <- function(theta, dat1, centerpersons)
{
    if (centerpersons){
        theta <- theta - stats::weighted.mean( theta, dat1[,2] )
    }
    return(theta)
}
