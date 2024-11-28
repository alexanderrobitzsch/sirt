## File Name: locpolycor_est_thresh_grad_fun.R
## File Version: 0.141


locpolycor_est_thresh_grad_fun <- function(x, y, x1, w, x0, model, K, eps=1e-10)
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
    dpred1 <- stats::dnorm(pred1)
    dpred2 <- stats::dnorm(pred2)
    if (model=='lin'){
        lpred1 <- stats::dnorm(pred1)*(x1-x0)
        lpred2 <- stats::dnorm(pred2)*(x1-x0)
    }

    arg <- ifelse( arg < eps, eps, arg )
    NP <- length(x)
    grad <- rep(NA, NP)
    for (pp in 1L:NP){
        # P(Y=0)=Phi(tau1) - 0
        # P(Y=1)=P(Y<=1) - P(Y=0)=Phi(tau2) - Phi(tau1), ...
        arg1 <- dpred1
        arg2 <- dpred2
        uu <- pp
        if (model=='lin' & (pp>(NP/2) )){
            arg1 <- lpred1
            arg2 <- lpred2
            uu <- pp - NP/2
        }
        # d log(arg) / d par=d par / arg
        a1 <- ifelse(y==uu-1, arg1, 0)
        a2 <- ifelse(y==uu, -arg2, 0)
        ad <- a1+a2
        grad[pp] <- - sum(w*ad/arg)
    }
    return(grad)
}
