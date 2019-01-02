## File Name: sirt_pem_create_parameter_index.R
## File Version: 0.09


sirt_pem_create_parameter_index <- function( parmlist )
{
    NP <- length(parmlist)
    parm_index <- list()
    parmlist_names <- names(parmlist)
    parm_index_names <- rep("", NP)
    last_index <- 0
    for (pp in 1:NP){
        x <- parmlist[[pp]]
        x_pp <- list()
        parm_index_names[pp] <- parmlist_names[pp]
        x_dim <- sirt_pem_extract_dimension(x=x)
        NX <- prod(x_dim)
        x_index <- last_index + seq( 1, NX )
        x_pp$index <- x_index
        last_index <- max(x_index)
        x_pp$dim <- x_dim
        parm_index[[pp]] <- x_pp
    }
    names(parm_index) <- parm_index_names
    #-- more meta-information about parameters
    parm_index[["__n_parameters"]] <- NP
    parm_index[["__length_parameter_vector"]] <- last_index
    return(parm_index)
}
