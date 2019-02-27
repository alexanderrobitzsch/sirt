## File Name: sirt_trim_increment.R
## File Version: 0.02

sirt_trim_increment <- function(increment, max_increment)
{
    eps <- 1E-10
    ci <- ceiling( abs(increment) / ( max_increment + eps ) )
    change <- ifelse( abs(increment) > max_increment, increment/(2*ci), increment )
    return(change)
}
