## File Name: linking_haberman_itempars_convert_process_matrices.R
## File Version: 0.01


linking_haberman_itempars_convert_process_matrices <- function(mat, est_pars)
{
    mat[ ! est_pars ] <- NA
    return(mat)
}
