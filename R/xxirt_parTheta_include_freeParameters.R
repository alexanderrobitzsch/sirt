## File Name: xxirt_parTheta_include_freeParameters.R
## File Version: 0.02


xxirt_parTheta_include_freeParameters <- function(customTheta, x)
{
    if (!is.null(x)){
        customTheta$par[ customTheta$est ] <- x
    }
    return(customTheta)
}
