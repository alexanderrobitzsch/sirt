## File Name: sirt_digamma1.R
## File Version: 0.01
## File Last Change: 2019-05-18


#-- derivative of digamma function
sirt_digamma1 <- function(x, h=1e-3)
{
    ( digamma(x+h) - digamma(x-h) ) / (2*h)
}
