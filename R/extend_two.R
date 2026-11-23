## File Name: extend_two.R
## File Version: 0.01

extend_two <- function(x)
{
    if (length(x)==1){
        x <- rep(x,2)
    }
    return(x)
}
