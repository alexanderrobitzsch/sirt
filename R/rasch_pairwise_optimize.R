## File Name: rasch_pairwise_optimize.R
## File Version: 0.222
## File Last Change: 2021-04-22


rasch_pairwise_optimize <- function(n.ij, n.ji, beta, zerosum, optimizer="nlminb",
        h=1e-4, maxiter=1000, estimator="MINCHI")
{
    I <- nrow(n.ij)
    ctl <- list(iter.max=maxiter)
    opt_fun1 <- function(x)
    {
        eps <- rasch_pairwise_compute_eps(x=x)
        t3 <- rasch_pairwise_optimize_opt_fun_terms(eps_horiz=eps,
                    eps_vert=eps, n.ij=n.ij, n.ji=n.ji)
        val <- sum(t3)
        return(val)
    }

    tol <- 1e-7
    if (estimator=="MINCHI"){
        y.ij <- n.ij / sqrt(n.ij+n.ji+tol)
        y.ji <- n.ji / sqrt(n.ij+n.ji+tol)
    }
    if (estimator=="ULS"){
        y.ij <- n.ij
        y.ji <- n.ji
    }
    opt_fun2 <- function(x)
    {
        eps <- rasch_pairwise_compute_eps(x=x)
        t3 <- rasch_pairwise_optimize_opt_fun_terms2(eps_horiz=eps,
                    eps_vert=eps, y.ij=y.ij, y.ji=y.ji, estimator=estimator)
        val <- sum(t3)
        return(val)
    }

    par_init <- beta - beta[1]
    par_init <- par_init[-1]

    grad_fun2 <- function(x)
    {
        eps <- rasch_pairwise_compute_eps(x=x)
        eps_h <- rasch_pairwise_compute_eps(x=x+h)
        eps_h2 <- rasch_pairwise_compute_eps(x=x-h)

        t2 <- rasch_pairwise_optimize_opt_fun_terms2(eps_horiz=eps_h,
                    eps_vert=eps, y.ij=y.ij, y.ji=y.ji, estimator=estimator)
        t3 <- rasch_pairwise_optimize_opt_fun_terms2(eps_horiz=eps,
                    eps_vert=eps_h, y.ij=y.ij, y.ji=y.ji, estimator=estimator)
        t21 <- rasch_pairwise_optimize_opt_fun_terms2(eps_horiz=eps_h2,
                    eps_vert=eps, y.ij=y.ij, y.ji=y.ji, estimator=estimator)
        t31 <- rasch_pairwise_optimize_opt_fun_terms2(eps_horiz=eps,
                    eps_vert=eps_h2, y.ij=y.ij, y.ji=y.ji, estimator=estimator)

        val2a <- rowSums(t2)
        val3a <- colSums(t3)
        val21 <- rowSums(t21)
        val31 <- colSums(t31)

        val2 <- val2a[-1]+val3a[-1]
        val1 <- val21[-1]+val31[-1]
        grad <- (val2-val1)/(2*h)
        return(grad)
    }

    res_opt <- sirt_optimizer(optimizer=optimizer, par=par_init, fn=opt_fun2,
                    grad=grad_fun2, control=ctl, hessian=FALSE)
    b <- c(0,-res_opt$par)
    if (zerosum){
        b <- b - mean(b)
    }
    eps <- exp(b)
    iter <- res_opt$iterations

    #--- output
    res <- list(eps=eps, b=b, iter=iter, res_opt=res_opt)
    return(res)
}
