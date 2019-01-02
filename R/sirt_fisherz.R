## File Name: sirt_fisherz.R
## File Version: 0.06


### just a copy of psych::fisherz
sirt_fisherz <- function(rho)
{
    0.5 * log((1 + rho)/(1 - rho))
}
