## File Name: rm_sdt_create_parm_index_modify_elements.R
## File Version: 0.10
## File Last Change: 2018-12-30


rm_sdt_create_parm_index_modify_elements <- function(x, start_index, type)
{
    is_fixed <- sum( x < 0 )
    ind2 <- start_index + match( x, unique(x[x>0]) ) - 1
    if (is.matrix(x) ){
        ind2 <- matrix(ind2, nrow=nrow(x), ncol=ncol(x))
    }
    x <- ind2
    if (is_fixed){
        x[ is.na(x) ] <- -9
    }
    if (is.matrix(x)){
        partable0 <- NULL
        x_dim <- dim(x)
        K <- x_dim[2]
        for (kk in 1:K){
            par1 <- data.frame(type=type, parindex=x[,kk], row=1:x_dim[1], col=kk)
            partable0 <- rbind(partable0, par1)
        }
    }
    if (is.vector(x)){
        x_dim <- dim(x)
        partable0 <- data.frame(type=type, parindex=x, row=seq_len(length(x)), col=NA)
    }
    partable0$est <- partable0$parindex > 0
    partable0$est <- partable0$est & ( ! duplicated(partable0$parindex) )
    partable0$fixed <- partable0$parindex < 0
    partable0$value <- NA
    return(partable0)
}
