## File Name: sirt_symmetrize.R
## File Version: 0.01
## File Last Change: 2022-01-27

sirt_symmetrize <- function(x)
{
    ( x+t(x) ) / 2
}
