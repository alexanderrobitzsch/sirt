## File Name: sirt_logdet.R
## File Version: 0.02

sirt_logdet <- function(x)
{
    as.numeric(determinant(x=x, logarithm=TRUE)$modulus)
}
