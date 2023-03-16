## File Name: sirt_symmetrize.R
## File Version: 0.01

sirt_symmetrize <- function(x)
{
    ( x+t(x) ) / 2
}
