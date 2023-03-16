## File Name: sirt_is_data.R
## File Version: 0.01
## File Last Change: 2022-01-12

sirt_is_data <- function(dat)
{
    is_data <- is.data.frame(dat) | is.matrix(dat)
    return(is_data)
}
