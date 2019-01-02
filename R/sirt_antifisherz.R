## File Name: sirt_antifisherz.R
## File Version: 0.03

## see psych::fisherz2r
sirt_antifisherz <- function(z)
{
    ( exp(2*z) - 1 ) / ( exp(2*z) + 1 )
}
