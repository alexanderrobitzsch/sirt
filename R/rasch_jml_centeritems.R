## File Name: rasch_jml_centeritems.R
## File Version: 0.04

rasch_jml_centeritems <- function(b, centeritems)
{
    if (centeritems){
        b <- b - mean(b)
    }
    return(b)
}
