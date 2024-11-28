## File Name: locpolycor_est_polycor_opt_fun.R
## File Version: 0.084


locpolycor_est_polycor_opt_fun <- function(x, w, x1, x0,
        thresh_ii1, thresh_ii2, thresh_jj1, thresh_jj2,
        model="const", eps=1e-10,
        package="pbivnorm", maxcor=0.99)
{
    requireNamespace(package)
    if (package=='pbivnorm'){
        pbvnorm_fun <- pbivnorm::pbivnorm
    }
    if (package=='pbv'){
        pbvnorm_fun <- pbv::pbvnorm
    }
    if (model=='const'){
        rho1 <- x
    }
    if (model=='lin'){
        rho1 <- x[1] + x[2] * (x1-x0)
        rho1 <- pmin( maxcor, pmax(-maxcor, rho1) )
    }

    d11 <- pbvnorm_fun(x=thresh_ii1, y=thresh_jj1, rho=rho1)
    d10 <- pbvnorm_fun(x=thresh_ii1, y=thresh_jj2, rho=rho1)
    d01 <- pbvnorm_fun(x=thresh_ii2, y=thresh_jj1, rho=rho1)
    d00 <- pbvnorm_fun(x=thresh_ii2, y=thresh_jj2, rho=rho1)
    val <- d11 - d10 - d01 + d00
    val <- ifelse(val<eps, eps, val)
    val <- - sum( w*log(val+eps) )
    return(val)
}
