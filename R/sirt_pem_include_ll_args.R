## File Name: sirt_pem_include_ll_args.R
## File Version: 0.01
## File Last Change: 2017-10-03 12:13:12

sirt_pem_include_ll_args <- function(ll_args, pem_parm, pem_pars, pem_parameter_index)
{
	for (pp in pem_pars){
		ll_args[[ pp ]] <- sirt_pem_extract_parameters( parm=pem_parm, parmgroup=pp, pem_parameter_index=pem_parameter_index )
	}								
	return(ll_args)
}
