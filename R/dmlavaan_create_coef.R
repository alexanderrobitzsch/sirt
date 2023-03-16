## File Name: dmlavaan_create_coef.R
## File Version: 0.01

dmlavaan_create_coef <- function(partable)
{
    coef1 <- partable$est1
    names(coef1) <- paste0(partable$parname, '_mod1')
    coef2 <- partable$est2
    names(coef2) <- paste0(partable$parname, '_mod2')
    coef <- c( coef1, coef2 )
    coef <- coef[ ! is.na(coef) ]
    return(coef)
}
