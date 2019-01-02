## File Name: noharm_sirt_outer_coefs.R
## File Version: 0.03

noharm_sirt_outer_coefs <- function(x)
{
    y <- as.matrix( TAM::tam_outer(x,x) )
    return(y)
}
