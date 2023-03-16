## File Name: sirt_fisherz.R
## File Version: 0.07
## File Last Change: 2019-04-30


### just a copy of fisherz from psych package
sirt_fisherz <- function(rho)
{
    0.5 * log((1 + rho)/(1 - rho))
}
