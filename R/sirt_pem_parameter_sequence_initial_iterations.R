## File Name: sirt_pem_parameter_sequence_initial_iterations.R
## File Version: 0.01
## File Last Change: 2017-10-03 11:37:14

sirt_pem_parameter_sequence_initial_iterations <- function( pem_parm, pem_parameter_sequence, iter )
{
	if (iter < 3){
		for (ii in 0:2){
			if (iter==ii){			
				pem_parameter_sequence[[ paste0("P",ii) ]] <- pem_parm
			}
		}
	}
	return(pem_parameter_sequence)
}
