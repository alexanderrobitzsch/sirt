## File Name: sirt_antifisherz.R
## File Version: 0.04

## see fisherz2r in psych package
sirt_antifisherz <- function(z)
{
    ( exp(2*z) - 1 ) / ( exp(2*z) + 1 )
}
