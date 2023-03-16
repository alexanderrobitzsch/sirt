## File Name: sirt_fisherz.R
## File Version: 0.07


### just a copy of fisherz from psych package
sirt_fisherz <- function(rho)
{
    0.5 * log((1 + rho)/(1 - rho))
}
