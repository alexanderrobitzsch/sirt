## File Name: noharm_sirt_preproc_parameter_table_matrix.R
## File Version: 0.05

noharm_sirt_preproc_parameter_table_matrix <- function(pattmat, valmat, patt_id,
    patt_label, minval, symm=FALSE)
{
    pattmat <- noharm_sirt_preproc_pattern_matrix(pattmat=pattmat, minval=minval,
                    symm=symm)
    I <- nrow(pattmat)
    D <- ncol(pattmat)
    parm1 <- data.frame(mat=patt_label, matid=patt_id, row=rep(1:I,D), col=rep(1:D, each=I) )
    parm1$index <- as.vector(pattmat)
    parm1$starts <- as.vector(valmat)
    parm1$est_par <- 1 - is.na(parm1$index)
    parm1$nonnull_par <- 1 - (1-parm1$est)*(parm1$starts==0)
    parm1$fixed <- 1*is.na(parm1$index)
    if (symm){
        parm1$est_par <- parm1$est_par * (parm1$row >=parm1$col )
    }
    return(parm1)
}
