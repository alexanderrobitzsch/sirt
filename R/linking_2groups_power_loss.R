## File Name: linking_2groups_power_loss.R
## File Version: 0.05

linking_2groups_power_loss <- function(x, pow, eps, deriv=FALSE)
{
    #--- no derivative
    if (!deriv){
        if (pow>0 & pow < 2){
            res <- ( x^2 + eps )^(pow/2)
        }
        if (pow==0){
            res <- x^2 / ( x^2 + eps )
        }
        if (pow==2){
            res <- x^2
        }
    } else {
    #--- derivative
        if (pow>0 & pow < 2){
            res <- x*pow*( x^2 + eps )^(pow/2-1)
        }
        if (pow==0){
            # res <- x^2 / ( x^2 + eps )
            res <- 2*x*eps  / (x^2 + eps)^2
        }
        if (pow==2){
            res <- 2*x
        }
    }
    return(res)
}
