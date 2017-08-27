## File Name: sirt_fisherz.R
## File Version: 0.01
## File Last Change: 2017-07-11 15:47:51


### just a copy of psych::fisherz
sirt_fisherz <- function(rho) 
{
    0.5 * log((1 + rho)/(1 - rho))
}
