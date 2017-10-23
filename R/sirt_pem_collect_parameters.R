## File Name: sirt_pem_collect_parameters.R
## File Version: 0.03
## File Last Change: 2017-10-03 11:15:33

sirt_pem_collect_parameters <- function( parmlist, pem_parameter_index )
{
	NV <- pem_parameter_index[["__length_parameter_vector"]]
	parm <- rep(0, NV)
	NP <- pem_parameter_index[["__n_parameters"]]
	NPL <- length(parmlist)
	for (pp in 1:NPL){
		var_pp <- names(parmlist)[pp]
		parm_index_pp <- pem_parameter_index[[ var_pp ]]
		x <- as.vector(parmlist[[pp]])
		parm[ parm_index_pp$index ] <- x
	}
	return(parm)
}
