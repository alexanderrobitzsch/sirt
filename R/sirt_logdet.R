## File Name: sirt_logdet.R
## File Version: 0.02
## File Last Change: 2022-01-25

sirt_logdet <- function(x)
{
    as.numeric(determinant(x=x, logarithm=TRUE)$modulus)
}
