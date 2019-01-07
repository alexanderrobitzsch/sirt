## File Name: noharm_sirt_optim_function_R.R
## File Version: 0.02


noharm_sirt_optim_function_R <- function(gamma_val, delta, I, wgtm, pm,
    b0.jk, b1.jk, b2.jk, b3.jk)
{
    val <- 0
    for (ii in 1:(I-1)){
        for (jj in (ii+1):I){
            if (wgtm[ii,jj] >0 ){
                x_ij <- gamma_val[ii,jj] / sqrt( delta[ii] * delta[jj] )
                pm_exp <- b0.jk[ii,jj] + b1.jk[ii,jj]*x_ij + b2.jk[ii,jj]*x_ij^2 + b3.jk[ii,jj]*x_ij^3
                val <- val + wgtm[ii,jj]*(pm[ii,jj] - pm_exp)^2
            }
        }
    }
    #-- output
    return(val)
}
