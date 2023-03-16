## File Name: sirt_pem_adjust_dimension.R
## File Version: 0.06
## File Last Change: 2018-12-30

sirt_pem_adjust_dimension <- function(x, x_dim )
{
    if ( length(x_dim)==2 ){
        x <- matrix(x, nrow=x_dim[1], ncol=x_dim[2] )
    }
    if ( length(x_dim) >=3){
        x <- array(x, dim=x_dim )
    }
    return(x)
}
