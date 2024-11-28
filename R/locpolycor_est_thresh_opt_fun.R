## File Name: locpolycor_est_thresh_opt_fun.R
## File Version: 0.153



locpolycor_est_thresh_opt_fun <- function(x, y, x1, w, x0, model, K, eps=1e-10)
{
    if (model=='const'){
        par1 <- c(-Inf,x,Inf)
        pred1 <- par1[y+2]
        pred2 <- par1[y+1]
    }
    if (model=='lin'){
        par1 <- c(-Inf,x[1L:K],Inf)
        par2 <- c(-Inf,x[K+(1L:K)],Inf)
        pred1 <- par1[y+2] + par2[y+2]*(x1-x0)
        pred1 <- ifelse(y==K, par1[y+2], pred1)
        pred2 <- par1[y+1] + par2[y+1]*(x1-x0)
        pred2 <- ifelse(y==0, par2[y+1], pred2)
    }
    arg <- stats::pnorm(pred1) - stats::pnorm(pred2)
    arg <- ifelse( arg < eps, eps, arg )
    val <- - sum(w*log(arg))
    return(val)
}
