## File Name: mgsem_cda_opt_evaluate_penalties.R
## File Version: 0.02

mgsem_cda_opt_evaluate_penalties <- function(x, opt_fun_args)
{
    res <- mgsem_evaluate_penalties(x=x,
                    partable=opt_fun_args$partable,
                    prior_list=opt_fun_args$prior_list,
                    technical=opt_fun_args$technical,
                    h=opt_fun_args$technical$h,
                    p=opt_fun_args$p_pen,
                    eps_approx=opt_fun_args$eps_approx,
                    deriv=FALSE, difflp_info=opt_fun_args$difflp_info,
                    loop_parms=opt_fun_args$loop_parms,
                    pen_type=opt_fun_args$pen_type, a_scad=opt_fun_args$a_scad)
    return(res)
}
