## File Name: xxirt_parTheta_extract_freeParameters.R
## File Version: 0.08
## File Last Change: 2018-12-30


xxirt_parTheta_extract_freeParameters <- function( customTheta )
{
    ind <- customTheta$est
    p1 <- customTheta$par[ ind ]
    names(p1) <- names(customTheta$par)[ind]
    return(p1)
}
