## File Name: sirt_display_function.R
## File Version: 0.02
## File Last Change: 2023-03-08

sirt_display_function <- function(length=66)
{
    disp <- paste0( rep('-',length), collapse='')
    cat(disp, '\n')
}
