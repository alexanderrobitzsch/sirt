## File Name: mgsem_add_increment.R
## File Version: 0.04
## File Last Change: 2022-02-25


mgsem_add_increment <- function(x, h, i1, i2=NULL, symm=FALSE )
{
    x1 <- x
    if (is.vector(x)){
        x1[i1] <- x[i1] + h
    } else {
        x1[i1,i2] <- x[i1,i2] + h
        if (symm & (i1!=i2) & (!is.null(i2)) ){
            x1[i2,i1] <- x1[i2,i1] + h
        }
    }
    return(x1)
}
