## File Name: plot.lsdm.R
## File Version: 0.02


plot.lsdm <- function(x, ...)
{
    graphics::matplot( x=x$theta, y=t(x$attr.curves), xlab=expression(theta),
        ylab="Attribute response curve", ylim=c(0,1), ... )
}
